"""T24 with the Chabrier conversion (see harness_T24.py). Usage: python3 harness_T24C.py <mask> <out> MOGFI@0.1"""
import sys
import harness_T24 as I
I.apply('C')
H, Y = I.H, I.Y
if __name__ == '__main__':
    Y.X.T.Cc.run(sys.argv[1], sys.argv[2], sys.argv[3:])
