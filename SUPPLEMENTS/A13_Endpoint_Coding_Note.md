# Supplementary Coding Calculations for A13

This optional note supplies the endpoint calculations behind A13's brief numerical illustration. It is separate from the main paper so the temporal argument does not depend on a lengthy discussion of clocks or encoding formats.

## Cosmological scale and fixed-width addresses

The age of our universe gives the choice a concrete scale. Taking 13.8 billion years as the interval, with 365.25 days per year, gives about $4.35\times10^{17}$ seconds.

| Resolution | Approximate number of positions | Fixed-width address |
|---|---:|---:|
| One second | $4.35\times10^{17}$ | 59 bits |
| One Planck time | $8.08\times10^{60}$ | 203 bits |

The second row uses $t_P\approx5.39\times10^{-44}$ seconds (NIST, 2022 CODATA). Planck time is the time light takes to travel one Planck length; using it here does not establish that physical time advances in such ticks. If 203 bits are additional to an otherwise unchanged description, the weight factor is $2^{-203}$, approximately $7.78\times10^{-62}$. The reduction is about $1.29\times10^{61}$, not $10^{203}$.

### What compression can change

A fixed-width address can select any position in the interval. A description of one particular position may be much shorter. For example, a rule can specify $2^{203}$ steps by encoding the base, exponent and operation. That expression does not specify our exact age; it demonstrates why a large integer need not have a comparably long shortest description. The measured age is itself approximate, rather than an exact Planck-time address.

The proposed two five-bit number fields do not settle the minimum. Even if each field represents values from 1 through 32, $32^{32}=2^{160}$ falls below the illustrated scale. This excludes that particular single-power format. It does not exclude short encodings using other operations or information already supplied by the model. An exact bit count requires a defined language. The next section gives a 16-bit stacked-exponent example: $2^{(3^5)}=2^{243}$. Sixteen bits is the minimum for reaching this scale within that stated grammar. This illustrates a duration cost; if all sixteen bits are additional to the reverse account, their weight factor is $2^{-16}=1/65{,}536$.

The same principle applies to an endpoint defined by an event. A rule such as “stop at the first observer satisfying criterion C” must specify the criterion and how it is evaluated. If that rule is effective and the event occurs, it can identify a time without listing a timestamp. Whether it costs more or less than a numerical address is a question for the complete encoding. Our existence does not establish that we are the first observers anywhere in the universe.

Most possible endpoints cannot all receive very short descriptions. For a fixed background, fewer than $2^m$ endpoints can have descriptions shorter than m bits, because there are only $2^m-1$ binary strings of those lengths. Among N uniformly sampled endpoints, the fraction with such descriptions is therefore less than $2^m/N$. This counting result does not establish the length of our particular endpoint, whose selection is tied to observer evidence.

## A sixteen-bit stacked-exponent illustration

For an explicit compressed integer language, let tag 0 introduce a positive integer literal in Elias gamma code and tag 1 introduce a power followed by two recursively encoded expressions. A literal n then uses $2\lfloor\log_2 n\rfloor+2$ bits including its tag. The expression $2^{203}$ uses $1+4+16=21$ bits, and $10^{61}$ uses $1+8+12=21$ bits. The shared decoder defines the language. These are exact lengths in this grammar, not universal minima or specifications of our exact endpoint.

Stacking improves the threshold example. Encode $2^{(3^5)}=2^{243}$ with two one-bit power tags and literal codes of four, four and six bits for 2, 3 and 5. The total is $1+4+1+4+6=16$ bits; the actual codeword is `1001010011000101`. This exceeds $2^{203}$ ticks. It specifies a duration beyond the cosmological threshold, not our exact measured age or the complete reversal procedure. Exponentiation is supplied by the shared decoder, so its implementation is not charged afresh at each call.

Sixteen is minimal for reaching that threshold in this grammar. Let $M(c)$ be the greatest integer expressible in exactly c bits. Its possibilities are a literal of that cost or $M(a)^{M(b)}$ with $a+b+1=c$. Exponentiation is monotone for positive integer arguments, so this recurrence covers every expression tree. Across all costs below sixteen, its largest value is $3^{31}=617{,}673{,}396{,}283{,}947$, far below the threshold. The sixteen-bit expression above reaches it. This supplies a checked lower bound within the declared format.

The bound is grammar-relative. Another prefix-free universal decoder can reserve `0` followed by eight bits for $2^n$, and `1p` for simulation of a fixed universal prefix-free machine on p. Its nine-bit code `011001011` produces $2^{203}$. Both functions and decoder are fixed before comparison. Thus sixteen is a shortest encoding in the worked grammar, not a minimum over all universal computers. The physical comparison must use one justified common language and count only the information the reverse account additionally needs.

## Reproducing the grammar bound

The recurrence check and its output are in the portfolio's EDITORIAL / ENDPOINT_CODING directory. Values are capped at the target without changing whether it is reached. The chosen clock and horizon remain illustrative; the reversal comparison concerns additional information.

## References

Grünwald, P. D., and Vitányi, P. M. B. (2008). Algorithmic information theory. https://arxiv.org/abs/0809.2754

National Institute of Standards and Technology. 2022 CODATA recommended values: Planck time. https://physics.nist.gov/cuu/Constants/Table/allascii.txt
