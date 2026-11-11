#!/usr/bin/env python3
"""TOE@home coordinator: hands out blocks of rules, checks every result against a second volunteer's, keeps the record.

Standard library only (Python 3.10 or later); the data is one SQLite file and a folder of result files.

  python3 toehome_server.py init   --db toe.sqlite [--classes b2_3,b2_4,b2_5,b3_4,b3_5,t2_3,b4_5]   (t2_4 can be added later)
  python3 toehome_server.py serve  --db toe.sqlite [--host 127.0.0.1] [--port 8742] [--quorum 2] [--behind-proxy] [--same-group-ok]
  python3 toehome_server.py report --db toe.sqlite            progress, the closest rules, contributors
  python3 toehome_server.py export --db toe.sqlite OUT.jsonl  every checked rule record, one per line
  python3 toehome_server.py load   --db toe.sqlite FILE...    result files from volunteers who ran offline
  python3 toehome_server.py trust  --db toe.sqlite CLIENT_ID  a machine of your own: its results count without a second one

What it trusts. Nothing a volunteer sends is run, and a result is accepted only if (1) it comes with the token issued for
that block, (2) it lists exactly the rules the block contains, which the coordinator works out itself, (3) every field has
the expected type, and (4) a second volunteer from another network address returns the same result: whole numbers equal,
readings equal to 0.002. A third volunteer breaks a disagreement. A rule that reads near four gets a long run by the same
two-volunteer procedure, and anything that passes the criterion there must still be rerun by hand before it is claimed.

What it keeps. A random label per volunteer machine, an optional nickname, and a salted hash of the address prefix (used
only so that the two results of a block come from different places). No address is stored.

Put it behind a web server with HTTPS (and pass --behind-proxy so the forwarded address is used). By default it listens on
127.0.0.1 only.
"""
import sys, os, json, time, hmac, hashlib, gzip, zlib, sqlite3, argparse, threading, secrets, html
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toehome as T

LEASE = 12 * 3600          # seconds a block stays reserved for a volunteer
MAX_OPEN = 64              # blocks one volunteer may hold at a time
MAX_BODY = 4 << 20         # largest upload, compressed
MAX_JSON = 32 << 20        # largest upload, unpacked
MAX_RESULTS = 5            # results per block before it is set aside for a look by hand
PER_MINUTE = 240           # requests per address group per minute

SCHEMA = '''
create table if not exists meta(key text primary key, value text);
create table if not exists units(uid text primary key, ord integer, state integer default 0, digest text, rules integer, done real);
create table if not exists leases(id integer primary key, uid text, client text, grp text, token text, issued real,
                                  returned real, digest text, verdict integer, file text, sec real);
create table if not exists clients(client text primary key, nick text, first real, last real, ok integer default 0,
                                   bad integer default 0, rules integer default 0, trusted integer default 0);
create table if not exists hits(uid text, i integer, rule text, label text, n integer, mm real, mid real, chain_d real,
                                z real, kind text, primary key(uid, i, label, kind));
create index if not exists leases_uid on leases(uid);
create index if not exists leases_client on leases(client, returned);
create index if not exists units_state on units(state, ord);
'''
LOCK = threading.Lock()

def connect(path):
    con = sqlite3.connect(path, timeout=60); con.row_factory = sqlite3.Row
    con.executescript(SCHEMA)
    return con

def meta(con, key, default=None):
    row = con.execute('select value from meta where key=?', (key,)).fetchone()
    return row['value'] if row else default

def cmd_init(a):
    con = connect(a.db)
    if meta(con, 'secret') is None:
        con.execute('insert into meta values(?,?)', ('secret', secrets.token_hex(32)))
        con.execute('insert into meta values(?,?)', ('salt', secrets.token_hex(16)))
        con.execute('insert into meta values(?,?)', ('version', str(T.VERSION)))
    if int(meta(con, 'version')) != T.VERSION: sys.exit('this database belongs to another version of the computation')
    n0 = con.execute('select coalesce(max(ord), 0) m from units').fetchone()['m']; n = n0
    for cls in a.classes.split(','):
        if cls not in T.CLASSES: sys.exit('unknown class ' + cls)
        for u in T.unit_list(cls):
            n += 1
            con.execute('insert or ignore into units(uid, ord) values(?,?)', (T.uid(u), n))
    con.commit()
    print('blocks in the database:', con.execute('select count(*) c from units').fetchone()['c'])

# ------------------------------------------------------------------ checking what a volunteer sends
def num(x, lo, hi, kind=(int,), none=False):
    if x is None and none: return None
    if isinstance(x, bool) or not isinstance(x, kind) or not lo <= x <= hi: raise ValueError('field out of range')
    return x

def interval(r):
    if r is None: return None
    if not isinstance(r, dict): raise ValueError('bad interval')
    out = {'diamonds': num(r.get('diamonds'), 1, 100), 'n': num(r.get('n'), 1, 10 ** 7), 'chain': num(r.get('chain'), 0, 10 ** 7, (int, float))}
    for key in ('f', 'mm', 'mid', 'chain_d', 'coarse_d'): out[key] = num(r.get(key), -1e6, 1e6, (int, float), none=True)
    return out

def sanitize(rec, long):
    """Rebuild a record from the fields this program knows, with their types checked. Anything else is dropped."""
    if not isinstance(rec, dict): raise ValueError('bad record')
    cap = T.LONG_EVENTS if long else T.EVENTS
    out = {'i': num(rec.get('i'), 0, 10 ** 12), 'rule': rec.get('rule'), 'ev': num(rec.get('ev'), 0, 40 * cap),
           'g': num(rec.get('g'), 0, T.LONG_GEN), 'rel': num(rec.get('rel'), 0, 10 ** 8), 'nodes': num(rec.get('nodes'), 0, 10 ** 9),
           'ops': num(rec.get('ops'), 0, 10 ** 10), 'cut': num(rec.get('cut'), 0, 1)}
    if not isinstance(out['rule'], str) or len(out['rule']) > 200: raise ValueError('bad rule text')
    fit = rec.get('fit') or {}
    if not isinstance(fit, dict): raise ValueError('bad fit')
    out['fit'] = {k: num(fit[k], -1e6, 1e6, (int, float), none=True) for k in ('exp_rate', 'power') if k in fit}
    if 'depth' in rec:
        out['depth'] = num(rec['depth'], 0, 10 ** 7); out['past_max'] = num(rec.get('past_max'), 0, 10 ** 8)
        out['max_causes'] = num(rec.get('max_causes'), 0, 10 ** 6); out['max_effects'] = num(rec.get('max_effects'), 0, 10 ** 8)
        if not isinstance(rec.get('tree'), bool): raise ValueError('bad field')
        out['tree'] = rec['tree']
        if 'cone_dim' in rec: out['cone_dim'] = num(rec['cone_dim'], -1e6, 1e6, (int, float), none=True)
        for lab in ('root', 'mid'):
            if rec.get(lab) is not None: out[lab] = interval(rec[lab])
        if long:
            ser = rec.get('series', [])
            if not isinstance(ser, list) or len(ser) > 16: raise ValueError('bad series')
            out['series'] = []
            for s in ser:
                if not isinstance(s, dict): raise ValueError('bad series')
                row = {'n': num(s.get('n'), 1, 10 ** 7), 'chain': num(s.get('chain'), 0, 10 ** 7), 'height': num(s.get('height'), 0, 10 ** 7)}
                for key in ('frac', 'f', 'mm', 'mid', 'chain_d', 'coarse_d'): row[key] = num(s.get(key), -1e6, 1e6, (int, float), none=True)
                out['series'].append(row)
    return out

def check_result(result):
    """Returns (unit, cleaned result). Raises ValueError unless the result lists exactly the rules of its block."""
    if not isinstance(result, dict) or result.get('v') != T.VERSION: raise ValueError('wrong version')
    u = T.parse_uid(result.get('uid')); long = u[0] == 'long'
    recs = result.get('rules')
    if not isinstance(recs, list): raise ValueError('no rule list')
    want = T.rules_of(u)
    if len(recs) != len(want): raise ValueError('the block holds %d rules, the result %d' % (len(want), len(recs)))
    A = T.engine(); L = T.lefts(u[1] if long else u[0])[u[2] if long else u[1]][0]
    clean = []
    for rec, (i, R) in zip(recs, want):
        rec = sanitize(rec, long)
        if rec['i'] != i or rec['rule'] != A.rule_str(L, R): raise ValueError('rule %d of the block does not match' % i)
        clean.append(rec)
    return u, {'v': T.VERSION, 'uid': T.uid(u), 'rules': clean}

# ------------------------------------------------------------------ handing out and judging
class Coordinator:
    def __init__(self, a):
        self.a = a; self.dir = os.path.splitext(a.db)[0] + '_results'; os.makedirs(self.dir, exist_ok=True)
        con = connect(a.db)
        self.secret = meta(con, 'secret'); self.salt = meta(con, 'salt')
        if self.secret is None: sys.exit('run "init" first')
        if int(meta(con, 'version')) != T.VERSION: sys.exit('this database belongs to another version of the computation')
        con.close(); self.rate = {}

    def group(self, address):
        """A salted hash of the address prefix: enough to tell two places apart, not enough to recover an address."""
        parts = address.split('.') if '.' in address else address.split(':')
        return hashlib.sha256((self.salt + '|' + ':'.join(parts[:3])).encode()).hexdigest()[:16]

    def allowed(self, grp):
        now = time.time(); t0, n = self.rate.get(grp, (now, 0))
        if now - t0 > 60: t0, n = now, 0
        self.rate[grp] = (t0, n + 1)
        if len(self.rate) > 100000: self.rate.clear()
        return n < PER_MINUTE

    def token(self, uid, client, issued):
        return hmac.new(self.secret.encode(), ('%s|%s|%.3f' % (uid, client, issued)).encode(), hashlib.sha256).hexdigest()[:32]

    def touch(self, con, client, nick):
        now = time.time()
        con.execute('insert or ignore into clients(client, nick, first, last) values(?,?,?,?)', (client, nick, now, now))
        con.execute('update clients set last=?, nick=case when ?<>\'\' then ? else nick end where client=?', (now, nick, nick, client))

    def work(self, msg, grp):
        client, nick = ident(msg)
        with LOCK:
            con = connect(self.a.db)
            try:
                self.touch(con, client, nick); now = time.time(); live = now - LEASE
                held = con.execute('select count(*) c from leases where client=? and returned is null and issued>?', (client, live)).fetchone()['c']
                if held >= MAX_OPEN: con.commit(); return {'wait': 600, 'msg': 'this machine already holds %d blocks' % held}
                left = con.execute('select count(*) c from units where state in (0, 2)').fetchone()['c']
                if not left: con.commit(); return {'done': True}
                trusted = con.execute('select trusted from clients where client=?', (client,)).fetchone()['trusted']
                start = -10 ** 18                               # long runs sit below zero, so they go out first
                while True:
                    rows = con.execute('''select u.uid, u.state, u.ord,
                        (select count(*) from leases l where l.uid=u.uid and l.returned is not null) nret,
                        (select count(*) from leases l where l.uid=u.uid and l.returned is null and l.issued>?) npend
                        from units u where u.state in (0, 2) and u.ord>? order by u.ord limit 300''', (live, start)).fetchall()
                    if not rows: con.commit(); return {'wait': 300, 'msg': 'every open block is with a volunteer right now'}
                    for r in rows:
                        need = self.a.quorum if r['state'] == 0 else r['nret'] + 1
                        if trusted: need = max(need, r['nret'] + 1) if r['state'] == 2 else 1
                        if r['nret'] + r['npend'] >= need: continue
                        mine = con.execute('select count(*) c from leases where uid=? and (client=? or (grp=? and ?=0)) and (returned is not null or issued>?)',
                                           (r['uid'], client, grp, int(self.a.same_group_ok or bool(trusted)), live)).fetchone()['c']
                        if mine: continue
                        tok = self.token(r['uid'], client, now)
                        con.execute('insert into leases(uid, client, grp, token, issued) values(?,?,?,?,?)', (r['uid'], client, grp, tok, now))
                        con.commit()
                        return {'uid': r['uid'], 'token': tok, 'lease': LEASE}
                    start = rows[-1]['ord']
            finally:
                con.close()

    def result(self, msg, grp, offline=False):
        client, nick = ident(msg)
        try: u, res = check_result(msg.get('result'))            # the slow part, outside the lock
        except ValueError as e: return {'status': 'refused', 'msg': str(e)}
        dg = T.digest(res); uid = res['uid']
        offline = offline or not msg.get('token')                # a block run without a reservation: one result like any other
        with LOCK:
            con = connect(self.a.db)
            try:
                self.touch(con, client, nick)
                unit = con.execute('select * from units where uid=?', (uid,)).fetchone()
                if unit is None: return {'status': 'refused', 'msg': 'this block is not part of the project'}
                if offline:
                    if con.execute('select 1 from leases where uid=? and client=?', (uid, client)).fetchone(): return {'status': 'duplicate'}
                    cur = con.execute('insert into leases(uid, client, grp, token, issued) values(?,?,?,?,?)', (uid, client, grp, '', time.time()))
                    lease_id = cur.lastrowid
                else:
                    lease = con.execute('select * from leases where uid=? and client=? and token=?', (uid, client, str(msg.get('token')))).fetchone()
                    if lease is None: return {'status': 'refused', 'msg': 'no block was issued to this machine under that token'}
                    if lease['returned'] is not None: return {'status': 'duplicate'}
                    lease_id = lease['id']
                path = os.path.join(self.dir, '%s__%d.json.gz' % (uid.replace(':', '_'), lease_id))
                with gzip.open(path, 'wt') as fh: json.dump(res, fh, separators=(',', ':'))
                sec = msg.get('sec') if isinstance(msg.get('sec'), (int, float)) and not isinstance(msg.get('sec'), bool) else None
                con.execute('update leases set returned=?, digest=?, file=?, sec=? where id=?', (time.time(), dg, path, sec, lease_id))
                status = self.judge(con, uid) if unit['state'] in (0, 2) else 'late'
                con.commit()
                return {'status': status}
            finally:
                con.close()

    def judge(self, con, uid):
        rows = con.execute('''select l.*, c.trusted from leases l join clients c on c.client=l.client
                              where l.uid=? and l.returned is not null order by l.id''', (uid,)).fetchall()
        groups = []                                              # results that agree with each other
        for r in rows:
            res = json.load(gzip.open(r['file'], 'rt'))
            for g in groups:
                if g['digest'] == r['digest'] and T.floats_agree(g['res'], res): g['rows'].append(r); break
            else: groups.append({'digest': r['digest'], 'res': res, 'rows': [r]})
        for g in groups:
            places = {r['grp'] for r in g['rows']} if not self.a.same_group_ok else {r['client'] for r in g['rows']}
            if len(places) >= self.a.quorum or any(r['trusted'] for r in g['rows']):
                n = len(g['res']['rules'])
                con.execute('update units set state=1, digest=?, rules=?, done=? where uid=?', (g['digest'], n, time.time(), uid))
                good = {r['id'] for r in g['rows']}
                for r in rows:
                    ok = r['id'] in good
                    con.execute('update leases set verdict=? where id=?', (int(ok), r['id']))
                    con.execute('update clients set ok=ok+?, bad=bad+?, rules=rules+? where client=?', (int(ok), int(not ok), n if ok else 0, r['client']))
                self.record_hits(con, uid, g['res'])
                return 'validated'
        if len(rows) >= MAX_RESULTS: con.execute('update units set state=3 where uid=?', (uid,)); return 'conflict'
        if len(rows) >= self.a.quorum and len(groups) > 1: con.execute('update units set state=2 where uid=?', (uid,)); return 'conflict'
        return 'waiting'

    def record_hits(self, con, uid, res):
        u = T.parse_uid(uid); long = u[0] == 'long'
        for rec in res['rules']:
            for lab, n, v, z in T.near_four(rec):
                con.execute('insert or ignore into hits values(?,?,?,?,?,?,?,?,?,?)', (uid, rec['i'], rec['rule'], lab, n, v[0], v[1], v[2], round(z, 2), 'long window' if long else 'window'))
                if not long:                                     # second stage: the same rule, run long, ahead of everything else
                    luid = T.uid(('long', u[0], u[1], u[2], rec['i']))
                    lo = con.execute('select coalesce(min(ord), 0) m from units').fetchone()['m']
                    con.execute('insert or ignore into units(uid, ord) values(?,?)', (luid, min(lo, 0) - 1))
            if long and T.candidate(rec)[0]:
                con.execute('insert or ignore into hits values(?,?,?,?,?,?,?,?,?,?)', (uid, rec['i'], rec['rule'], 'series', 0, None, None, None, None, 'CANDIDATE'))

    def stats(self):
        con = connect(self.a.db)
        try:
            out = {'version': T.VERSION, 'classes': [], 'rules_checked': 0}
            for cls in ('long',) + T.ORDER:
                row = con.execute('''select count(*) n, coalesce(sum(state=1), 0) ok, coalesce(sum(state in (2, 3)), 0) bad, coalesce(sum(case when state=1 then rules end), 0) rules
                                     from units where uid like ?''', (cls + ':%',)).fetchone()
                if row['n']: out['classes'].append({'class': cls, 'blocks': row['n'], 'checked': row['ok'], 'in_dispute': row['bad'], 'rules_checked': row['rules']})
                if cls != 'long': out['rules_checked'] += row['rules']
            out['near_four'] = [dict(r) for r in con.execute('select rule, label, n, mm, mid, chain_d, z, kind from hits order by kind, z limit 40')]
            out['candidates'] = con.execute("select count(*) c from hits where kind='CANDIDATE'").fetchone()['c']
            out['volunteers'] = con.execute('select count(*) c from clients where ok>0').fetchone()['c']
            out['top'] = [{'who': r['nick'] or 'anonymous ' + r['client'][:6], 'blocks': r['ok'], 'rules': r['rules']}
                          for r in con.execute('select * from clients where ok>0 order by rules desc limit 25')]
            return out
        finally:
            con.close()

def ident(msg):
    client = msg.get('client'); nick = msg.get('nick') or ''
    if not isinstance(client, str) or not 8 <= len(client) <= 64 or not all(ch in '0123456789abcdef' for ch in client): raise ValueError('bad machine label')
    if msg.get('v') != T.VERSION: raise ValueError('this coordinator runs version %d of the computation; please download the current client' % T.VERSION)
    nick = ''.join(ch for ch in str(nick) if ch.isalnum() or ch in ' _-')[:24]
    return client, nick

def page(st):
    e = html.escape
    rows = ''.join('<tr><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td></tr>' % (e(c['class']), c['blocks'], c['checked'], c['in_dispute'], c['rules_checked']) for c in st['classes'])
    near = ''.join('<tr><td><code>%s</code></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>' % tuple(e(str(r[k])) for k in ('rule', 'kind', 'n', 'mm', 'mid', 'chain_d', 'z')) for r in st['near_four'])
    top = ''.join('<tr><td>%s</td><td>%d</td><td>%d</td></tr>' % (e(t['who']), t['blocks'], t['rules']) for t in st['top'])
    return ('<!doctype html><meta charset="utf-8"><title>TOE@home</title><style>body{font:16px/1.5 system-ui,sans-serif;max-width:60rem;margin:2rem auto;padding:0 1rem}'
            'table{border-collapse:collapse;margin:1rem 0}td,th{border:1px solid #bbb;padding:.3rem .6rem;text-align:left}</style>'
            '<h1>TOE@home</h1><p>A search for a short rule that grows four-dimensional spacetime. %d rules checked by %d volunteers. '
            'Rules passing the four-dimension criterion on a long run: <b>%d</b> (each must be rerun by hand before it is claimed).</p>'
            '<h2>Progress</h2><table><tr><th>Class</th><th>Blocks</th><th>Checked twice</th><th>In dispute</th><th>Rules checked</th></tr>%s</table>'
            '<h2>Rules reading between 3.5 and 4.5 on all three tests</h2><table><tr><th>Rule</th><th>Stage</th><th>Elements</th><th>Ordering</th><th>Midpoint</th><th>Chain</th><th>Distance from four (sd)</th></tr>%s</table>'
            '<h2>Contributors</h2><table><tr><th>Who</th><th>Blocks</th><th>Rules</th></tr>%s</table>') % (st['rules_checked'], st['volunteers'], st['candidates'], rows, near or '<tr><td colspan=7>none so far</td></tr>', top)

class Handler(BaseHTTPRequestHandler):
    server_version = 'toehome'; co = None
    def log_message(self, fmt, *args): pass                      # no addresses in any log
    def address(self):
        if self.co.a.behind_proxy:
            fwd = self.headers.get('X-Forwarded-For', '')
            if fwd: return fwd.split(',')[-1].strip()
        return self.client_address[0]
    def send(self, code, obj, kind='application/json'):
        body = obj.encode() if isinstance(obj, str) else json.dumps(obj).encode()
        self.send_response(code); self.send_header('Content-Type', kind + '; charset=utf-8')
        self.send_header('Content-Length', str(len(body))); self.send_header('X-Content-Type-Options', 'nosniff'); self.end_headers()
        self.wfile.write(body)
    def do_GET(self):
        grp = self.co.group(self.address())
        if not self.co.allowed(grp): return self.send(429, {'msg': 'too many requests'})
        if self.path == '/stats.json': return self.send(200, self.co.stats())
        if self.path in ('/', '/index.html'): return self.send(200, page(self.co.stats()), 'text/html')
        self.send(404, {'msg': 'not found'})
    def do_POST(self):
        grp = self.co.group(self.address())
        if not self.co.allowed(grp): return self.send(429, {'msg': 'too many requests'})
        try:
            n = int(self.headers.get('Content-Length', '0'))
            if not 0 < n <= MAX_BODY: return self.send(413, {'msg': 'upload too large'})
            raw = self.rfile.read(n)
            if self.headers.get('Content-Encoding') == 'gzip':
                d = zlib.decompressobj(31); raw = d.decompress(raw, MAX_JSON)
                if d.unconsumed_tail: return self.send(413, {'msg': 'upload too large'})
            msg = json.loads(raw.decode())
            if not isinstance(msg, dict): raise ValueError('bad request')
            if self.path == '/work': return self.send(200, self.co.work(msg, grp))
            if self.path == '/result': return self.send(200, self.co.result(msg, grp))
            self.send(404, {'msg': 'not found'})
        except (ValueError, zlib.error, UnicodeDecodeError) as e:
            self.send(400, {'msg': str(e)[:200]})
        except Exception:
            self.send(500, {'msg': 'server error'})

def cmd_serve(a):
    Handler.co = Coordinator(a)
    srv = ThreadingHTTPServer((a.host, a.port), Handler); srv.daemon_threads = True
    print('TOE@home coordinator, version %d, on http://%s:%d (stop with Ctrl-C)' % (T.VERSION, a.host, a.port), flush=True)
    try: srv.serve_forever()
    except KeyboardInterrupt: pass

def cmd_report(a):
    a.quorum = 2; a.same_group_ok = False; a.behind_proxy = False
    st = Coordinator(a).stats()
    print('rules checked: %d   volunteers: %d   passing the criterion on a long run: %d' % (st['rules_checked'], st['volunteers'], st['candidates']))
    for c in st['classes']: print('  %-5s blocks %6d   checked %6d   in dispute %4d   rules %8d' % (c['class'], c['blocks'], c['checked'], c['in_dispute'], c['rules_checked']))
    print('rules reading 3.5 to 4.5 on all three tests:', len(st['near_four']) or 'none')
    for r in st['near_four']: print('  %-12s %s  n=%s  %s / %s / %s  distance %s sd' % (r['kind'], r['rule'], r['n'], r['mm'], r['mid'], r['chain_d'], r['z']))
    for t in st['top']: print('  %-30s blocks %5d  rules %7d' % (t['who'], t['blocks'], t['rules']))

def cmd_export(a):
    con = connect(a.db); n = 0
    with open(a.out, 'w') as fh:
        for u in con.execute('select uid, digest from units where state=1 order by ord'):
            row = con.execute('select file from leases where uid=? and digest=? and verdict=1 order by id limit 1', (u['uid'], u['digest'])).fetchone()
            for rec in json.load(gzip.open(row['file'], 'rt'))['rules']:
                fh.write(json.dumps(dict(rec, block=u['uid'])) + '\n'); n += 1
    print(n, 'checked rule records written to', a.out)

def cmd_load(a):
    a.quorum = getattr(a, 'quorum', 2); a.same_group_ok = False; a.behind_proxy = False
    co = Coordinator(a)
    for path in a.files:
        try:
            pack = json.load(gzip.open(path, 'rt'))
            msg = {'v': T.VERSION, 'client': pack.get('client'), 'nick': '', 'result': pack.get('result'), 'sec': pack.get('sec')}
            print(os.path.basename(path), co.result(msg, 'file:' + str(pack.get('client'))[:12], offline=True))
        except (OSError, ValueError, EOFError) as e:
            print(os.path.basename(path), 'not read:', e)

def cmd_trust(a):
    con = connect(a.db); now = time.time()
    con.execute('insert or ignore into clients(client, nick, first, last) values(?,?,?,?)', (a.client, '', now, now))
    con.execute('update clients set trusted=1 where client=?', (a.client,)); con.commit()
    print('trusted:', a.client)

def main():
    ap = argparse.ArgumentParser(description='TOE@home coordinator'); sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('init'); p.add_argument('--db', required=True); p.add_argument('--classes', default=','.join(T.ORDER[:-1])); p.set_defaults(f=cmd_init)
    p = sub.add_parser('serve'); p.add_argument('--db', required=True); p.add_argument('--host', default='127.0.0.1'); p.add_argument('--port', type=int, default=8742)
    p.add_argument('--quorum', type=int, default=2); p.add_argument('--behind-proxy', action='store_true'); p.add_argument('--same-group-ok', action='store_true'); p.set_defaults(f=cmd_serve)
    p = sub.add_parser('report'); p.add_argument('--db', required=True); p.set_defaults(f=cmd_report)
    p = sub.add_parser('export'); p.add_argument('--db', required=True); p.add_argument('out'); p.set_defaults(f=cmd_export)
    p = sub.add_parser('load'); p.add_argument('--db', required=True); p.add_argument('files', nargs='+'); p.set_defaults(f=cmd_load)
    p = sub.add_parser('trust'); p.add_argument('--db', required=True); p.add_argument('client'); p.set_defaults(f=cmd_trust)
    a = ap.parse_args(); a.f(a)

if __name__ == '__main__':
    main()
