Lab 5 - Numerical Interpolation (Newton Forward & Backward)
Run:   python3 interpolation.py
Test:  python3 interpolation.py < sample_input.txt

Input : n, x values (equally spaced), f(x) values, then query values (blank line to stop).
Output: difference table; for each query -> method (FORWARD/BACKWARD), origin x0, h, p,
        term-by-term calculation, interpolated value; finally the max-degree polynomial.
Rule  : first half of table -> forward, origin = point just before x
        second half         -> backward, origin = point just after x   (so |p| < 1)
Values may be integers, decimals or fractions like 1/3 (exact arithmetic is used).
