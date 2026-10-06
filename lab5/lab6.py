
from fractions import Fraction
from math import factorial, floor, ceil

LINE = "=" * 78
THIN = "-" * 78


# ------------------------------------------------------------------ helpers
def num(s):
    return Fraction(s.strip())            # accepts 2, 2.5, 1/3


def fmt(v, digits=7):
    v = Fraction(v)
    if v.denominator == 1:
        return str(v.numerator)
    return f"{float(v):.{digits}g}"


def read_list(label, n):
    vals = []
    for index in range(n):
        while True:
            try:
                vals.append(num(input(f"{label} {index + 1}/{n} : ")))
                break
            except (ValueError, ZeroDivisionError):
                print("  ! Invalid number - please enter this value again.")
    return vals


def equally_spaced(xs):
    if len(xs) < 2:
        return True
    h = xs[1] - xs[0]
    return h != 0 and all(xs[k + 1] - xs[k] == h for k in range(len(xs) - 1))


def clamp(i, lo, hi):
    return max(lo, min(i, hi))


def binom(a, k):
    r = Fraction(1)
    for j in range(k):
        r = r * (a - j) / (j + 1)
    return r


def sup(k):
    return "" if k == 1 else {2: "²", 3: "³"}.get(k, f"^{k}")


# ------------------------------------------------------------------ difference table
def build_table(f):
    """table[i][k] = Delta^k f_i  (forward difference of order k at index i)."""
    n = len(f)
    t = [[f[i]] for i in range(n)]
    for k in range(1, n):
        for i in range(n - k):
            t[i].append(t[i + 1][k - 1] - t[i][k - 1])
    return t


def last_order(xs, t):
    n = len(xs)
    return max(k for k in range(n)
               if k == 0 or any(len(r) > k and abs(float(r[k])) > 1e-12 for r in t))


def mark(s, used_flag):
    return f"[{s}]" if used_flag else s


def print_table(xs, t, used=None, title="DIFFERENCE TABLE (forward differences)"):
    """forward-difference table.  used = set of (r,k) = Delta^k f_r cells on the PATH (shown in [ ])"""
    used = used or set()
    n = len(xs)
    last = last_order(xs, t)
    head = ["x", "f(x)"] + ["Δ" + ("" if k == 1 else "²" if k == 2 else "³" if k == 3 else f"^{k}") + "f"
                               for k in range(1, last + 1)]
    rows = [[fmt(xs[i])] + [mark(fmt(v), (i, k) in used) for k, v in enumerate(t[i][:last + 1])]
            for i in range(n)]
    w = max(len(s) for r in rows + [head] for s in r) + 3
    print("\n" + title)
    print(THIN)
    print("".join(s.rjust(w) for s in head))
    print(THIN)
    for r in rows:
        print("".join(s.rjust(w) for s in r))
    print(THIN)
    if last < n - 1:
        print(f"(differences of order {last + 1} and higher are all zero)")


def print_central_table(xs, t, used=None,
                        title="CENTRAL DIFFERENCE TABLE  (δ^k f placed midway between the entries it is made from)"):
    """central layout: Delta^k f_r drawn at row 2r+k.  used cells are shown in [ ]"""
    used = used or set()
    n = len(xs)
    last = last_order(xs, t)
    head = ["x", "f(x)"] + ["δ" + sup(k) + "f" for k in range(1, last + 1)]
    rows = []
    for r in range(2 * n - 1):
        line = [""] * (last + 2)
        if r % 2 == 0:
            line[0] = fmt(xs[r // 2])
        for k in range(0, last + 1):
            if (r - k) % 2 == 0:
                i = (r - k) // 2
                if 0 <= i < n and k < len(t[i]):
                    line[k + 1] = mark(fmt(t[i][k]), (i, k) in used)
        rows.append(line)
    w = max(len(s) for r in rows + [head] for s in r) + 3
    print("\n" + title)
    print(THIN)
    print("".join(s.rjust(w) for s in head))
    print(THIN)
    for r in rows:
        print("".join(s.rjust(w) for s in r))
    print(THIN)
    if last < n - 1:
        print(f"(differences of order {last + 1} and higher are all zero)")


def cell_name(xs, h, t, r, k):
    """absolute name + value of Delta^k f_r  =  delta^k f at x = x_r + k h/2"""
    nm = "f" if k == 0 else f"δ{sup(k)}f"
    return f"{nm}(x={fmt(xs[r] + k * h / 2)})={fmt(t[r][k])}"


# ================================================================== NEWTON FORWARD / BACKWARD (original)
def select(xs, h, q, force=None):
    n = len(xs)
    i = floor((q - xs[0]) / h)            # tabular point just before q
    i = max(0, min(i, n - 1))
    fwd_terms = n - 1 - i
    bwd_terms = i
    if force:
        method = force
    else:
        method = "FORWARD" if fwd_terms >= bwd_terms else "BACKWARD"
    return method, i, fwd_terms, bwd_terms


def interpolate(xs, t, h, q, force=None):
    n = len(xs)
    method, i, fo, bo = select(xs, h, q, force)
    p = (q - xs[i]) / h
    if method == "FORWARD":
        diffs = [t[i][k] for k in range(n - i)]          # Delta^k f_i
    else:
        diffs = [t[i - k][k] for k in range(i + 1)]      # Nabla^k f_i = Delta^k f_(i-k)
    coefs = []
    c = Fraction(1)
    for k in range(len(diffs)):
        if k > 0:
            c = c * (p - (k - 1) if method == "FORWARD" else p + (k - 1)) / k
        coefs.append(c)
    total = sum(c * d for c, d in zip(coefs, diffs))
    return method, i, fo, bo, p, diffs, coefs, total


def term_text(k, p, d, method):
    """lecture style: (1/k!) p(p-1)...  x  difference"""
    if k == 0:
        return fmt(d)
    sgn = "-" if method == "FORWARD" else "+"
    fac = "".join(f"({fmt(p)}{sgn if j else ''}{j if j else ''})" if j else f"({fmt(p)})"
                  for j in range(k))
    pre = f"1/{k}! " if k > 1 else ""
    return f"{pre}{fac}({fmt(d)})"


def show_query(xs, t, h, q, force=None):
    n = len(xs)
    method, i, fo, bo, p, diffs, coefs, total = interpolate(xs, t, h, q, force)
    print("\n" + LINE)
    print(f"  QUERY  x = {fmt(q)}")
    print(LINE)
    lo = max(0, min(i, n - 2))
    print(f"  Table range          : [{fmt(xs[0])}, {fmt(xs[-1])}],  h = {fmt(h)}")
    print(f"  Query lies between   : {fmt(xs[lo])} and {fmt(xs[lo + 1])}")
    print(f"  Differences at origin: forward gives up to order {fo}, backward up to order {bo}")
    if force:
        better = "FORWARD" if fo > bo else "BACKWARD" if bo > fo else None
        why = "selected from the menu"
        if better and better != force:
            why += f"  (note: {better} would have more differences here)"
    else:
        why = ("more forward differences available" if fo > bo else
               "more backward differences available" if bo > fo else
               "equal number of differences -> forward chosen")
    print(f"  METHOD               : NEWTON {method}   ({why})")
    print(f"  ORIGIN               : x0 = {fmt(xs[i])}   (tabular point just before x)")
    print(f"  p = (x - x0) / h     : ({fmt(q)} - {fmt(xs[i])}) / {fmt(h)} = {fmt(p)}"
          f"   {'(0 <= p < 1  OK)' if 0 <= p < 1 else '(outside 0..1)'}")
    print(THIN)
    if method == "FORWARD":
        print("  FORMULA: f = f0 + p*Δf0 + p(p-1)/2! * Δ²f0 + p(p-1)(p-2)/3! * Δ³f0 + ...")
    else:
        print("  FORMULA: f = f0 + p*∇f0 + p(p+1)/2! * ∇²f0 + p(p+1)(p+2)/3! * ∇³f0 + ...")
    name = "Δ" if method == "FORWARD" else "∇"
    print(f"  values at origin x0 = {fmt(xs[i])}:  " +
          ",  ".join(f"{name}{k if k > 1 else ''}f0 = {fmt(d)}" if k else f"f0 = {fmt(d)}"
                      for k, d in enumerate(diffs)))
    print(THIN)
    cells = [(i, k) for k in range(len(diffs))] if method == "FORWARD" else [(i - k, k) for k in range(len(diffs))]
    print("  PATH FOLLOWED (" + ("straight DOWN the diagonal starting at x0: f0 -> Δf0 -> Δ²f0 -> ..."
                                  if method == "FORWARD" else
                                  "straight UP the diagonal ending at x0: f0 -> ∇f0 -> ∇²f0 -> ...") + ")")
    sym = "Δ" if method == "FORWARD" else "∇"
    print("   " + "  ->  ".join((f"f0={fmt(t[r][k])}" if k == 0 else f"{sym}{sup(k)}f0={fmt(t[r][k])}")
                                for r, k in cells) + f"      (f0 means f at x0 = {fmt(xs[i])})")
    print_table(xs, t, set(cells), "  Table with the path marked in [ ]  (forward-difference layout)")
    print("  SUBSTITUTION:")
    print("   f = " + "\n     + ".join(term_text(k, p, d, method) for k, d in enumerate(diffs)))
    print("\n  TERM BY TERM:")
    for k, (c, d) in enumerate(zip(coefs, diffs)):
        print(f"   term {k}:  {fmt(c):>12} x {fmt(d):>10}  =  {fmt(c * d):>12}")
    print("  " + "-" * 44)
    print(f"   TOTAL                              =  {fmt(total, 10):>12}")
    print(LINE)
    print(f"  RESULT:  f({fmt(q)}) = {fmt(total, 10)}   [method: NEWTON {method}, origin x0 = {fmt(xs[i])}]")
    print(LINE)


# ================================================================== CENTRAL FORMULAE
# each returns (name, origin_index, p, terms)   terms = [(label, coef, diff), ...]

def gauss_forward(xs, t, h, x):
    n = len(xs)
    i = clamp(floor((x - xs[0]) / h), 0, n - 1)
    p = (x - xs[i]) / h
    terms = [("f0", Fraction(1), t[i][0], [(i, 0)])]
    for k in range(1, n):
        r = i - k // 2 if k % 2 == 0 else i - (k - 1) // 2
        if r < 0 or r + k > n - 1:
            break
        lab = f"δ{sup(k)}f0" if k % 2 == 0 else f"δ{sup(k)}f(1/2)"
        terms.append((lab, binom(p + (k - 1) // 2, k), t[r][k], [(r, k)]))
    return "Gauss Forward", i, p, terms


def gauss_backward(xs, t, h, x):
    n = len(xs)
    i = clamp(ceil((x - xs[0]) / h), 0, n - 1)
    p = (x - xs[i]) / h
    terms = [("f0", Fraction(1), t[i][0], [(i, 0)])]
    for k in range(1, n):
        r = i - k // 2 if k % 2 == 0 else i - (k - 1) // 2 - 1
        if r < 0 or r + k > n - 1:
            break
        lab = f"δ{sup(k)}f0" if k % 2 == 0 else f"δ{sup(k)}f(-1/2)"
        terms.append((lab, binom(p + k // 2, k), t[r][k], [(r, k)]))
    return "Gauss Backward", i, p, terms


def stirling(xs, t, h, x):
    n = len(xs)
    i = clamp(floor((x - xs[0]) / h + Fraction(1, 2)), 0, n - 1)     # nearest point
    p = (x - xs[i]) / h
    p2 = p * p
    terms = [("f0", Fraction(1), t[i][0], [(i, 0)])]
    for k in range(1, n):
        if k % 2 == 1:
            m = (k - 1) // 2
            lo, hi = i - m - 1, i - m
            if lo < 0 or hi + k > n - 1:
                break
            c = p
            for j in range(1, m + 1):
                c *= (p2 - j * j)
            c /= factorial(k)
            terms.append((f"½[δ{sup(k)}f(-1/2)+δ{sup(k)}f(1/2)]", c, (t[lo][k] + t[hi][k]) / 2, [(lo, k), (hi, k)]))
        else:
            m = k // 2
            r = i - m
            if r < 0 or r + k > n - 1:
                break
            c = p2
            for j in range(1, m):
                c *= (p2 - j * j)
            c /= factorial(k)
            terms.append((f"δ{sup(k)}f0", c, t[r][k], [(r, k)]))
    return "Stirling", i, p, terms


def bessel(xs, t, h, x):
    n = len(xs)
    i = clamp(floor((x - xs[0]) / h), 0, max(n - 2, 0))
    p = (x - xs[i]) / h
    terms = [("f0", Fraction(1), t[i][0], [(i, 0)])]
    if n > 1:
        terms.append(("δf(1/2)", p, t[i][1], [(i, 1)]))
    for k in range(2, n):
        m = k // 2
        if m > i or i + m + 1 > n - 1:
            break
        if k % 2 == 0:
            terms.append((f"½[δ{sup(k)}f0+δ{sup(k)}f1]", binom(p + m - 1, k),
                          (t[i - m][k] + t[i - m + 1][k]) / 2, [(i - m, k), (i - m + 1, k)]))
        else:
            terms.append((f"δ{sup(k)}f(1/2)", (p - Fraction(1, 2)) * binom(p + m - 1, k - 1) / k,
                          t[i - m][k], [(i - m, k)]))
    return "Bessel", i, p, terms


def everett(xs, t, h, x):
    n = len(xs)
    i = clamp(floor((x - xs[0]) / h), 0, max(n - 2, 0))
    p = (x - xs[i]) / h
    qc = 1 - p
    terms = []
    m = 0
    while m <= i and i + m + 1 <= n - 1:
        l0 = "f0" if m == 0 else f"δ{sup(2 * m)}f0"
        l1 = "f1" if m == 0 else f"δ{sup(2 * m)}f1"
        terms.append((l0, binom(qc + m, 2 * m + 1), t[i - m][2 * m], [(i - m, 2 * m)]))
        terms.append((l1, binom(p + m, 2 * m + 1), t[i - m + 1][2 * m], [(i - m + 1, 2 * m)]))
        m += 1
    return "Everett", i, p, terms


# key -> (function, formula text, suited p-range text, p-in-range test)
CENTRAL = {
    "1": (stirling, "f = f0 + p/2[δf(-1/2)+δf(1/2)] + p²/2! δ²f0 + p(p²-1)/3! ·½[δ³f(-1/2)+δ³f(1/2)] + p²(p²-1)/4! δ⁴f0 + ...",
          "-0.25 < p < 0.25 (origin = nearest point)", lambda p: -Fraction(1, 4) < p < Fraction(1, 4)),
    "2": (bessel, "f = f0 + pδf(1/2) + p(p-1)/2! ·½[δ²f0+δ²f1] + p(p-1)(p-1/2)/3! δ³f(1/2) + (p+1)p(p-1)(p-2)/4! ·½[δ⁴f0+δ⁴f1] + ...",
          "p near 0.5 (0.25 .. 0.75)", lambda p: Fraction(1, 4) <= p <= Fraction(3, 4)),
    "3": (everett, "f = q f0 + q(q²-1)/3! δ²f0 + q(q²-1)(q²-4)/5! δ⁴f0 + ... + p f1 + p(p²-1)/3! δ²f1 + p(p²-1)(p²-4)/5! δ⁴f1 + ...   (q = 1-p)",
          "0 <= p < 1", lambda p: 0 <= p < 1),
    "4": (gauss_forward, "f = f0 + pδf(1/2) + p(p-1)/2! δ²f0 + (p+1)p(p-1)/3! δ³f(1/2) + (p+1)p(p-1)(p-2)/4! δ⁴f0 + ...",
          "0 <= p < 1", lambda p: 0 <= p < 1),
    "5": (gauss_backward, "f = f0 + pδf(-1/2) + (p+1)p/2! δ²f0 + (p+1)p(p-1)/3! δ³f(-1/2) + (p+2)(p+1)p(p-1)/4! δ⁴f0 + ...",
          "-1 < p <= 0 (origin = point just after x)", lambda p: -1 < p <= 0),
}


def max_order(terms):
    best = 0
    for lab, *_rest in terms:
        for k in range(2, 10):
            if "δ" + sup(k) in lab:
                best = max(best, k)
        if "δf" in lab:
            best = max(best, 1)
    return best


PATH_DESC = {
    "1": "zig-zag through the MIDDLE row x0:  f0 -> mean of the two δ (above & below x0) -> δ²f0 -> mean of the two δ³ -> δ⁴f0 -> ...",
    "2": "zig-zag BETWEEN the rows x0 and x1:  f0 -> δf(1/2) -> mean of δ²f0, δ²f1 -> δ³f(1/2) -> mean of δ⁴f0, δ⁴f1 -> ...",
    "3": "TWO straight horizontal lines of EVEN differences, one along row x0 (f0, δ²f0, δ⁴f0 ...) and one along row x1 (f1, δ²f1, δ⁴f1 ...)",
    "4": "zig-zag starting at x0 going DOWN first:  f0 -> δf(1/2) -> δ²f0 -> δ³f(1/2) -> δ⁴f0 -> ...",
    "5": "zig-zag starting at x0 going UP first:  f0 -> δf(-1/2) -> δ²f0 -> δ³f(-1/2) -> δ⁴f0 -> ...",
}


def show_path(xs, t, h, key, terms, show_table):
    print("  PATH FOLLOWED: " + PATH_DESC[key])
    if key == "3":
        for nm, ser in (("row x0", terms[0::2]), ("row x1", terms[1::2])):
            print(f"   {nm}: " + "  ->  ".join(cell_name(xs, h, t, *tm[3][0]) for tm in ser))
    else:
        parts = []
        for tm in terms:
            cs = [cell_name(xs, h, t, r, k) for r, k in tm[3]]
            parts.append(cs[0] if len(cs) == 1 else "mean{ " + " , ".join(cs) + " }")
        print("   " + "  ->  ".join(parts))
    if show_table:
        used = set(c for tm in terms for c in tm[3])
        print_central_table(xs, t, used, "  Table with the path marked in [ ]  (central layout)")


def show_central_detail(xs, t, h, x, key, res, show_table):
    name, i, p, terms, total, ftxt, rng, ok = res
    print("\n" + THIN)
    print(f"  {name.upper()}   origin x0 = {fmt(xs[i])},  "
          f"p = ({fmt(x)} - {fmt(xs[i])})/{fmt(h)} = {fmt(p)}")
    print(f"  suited to: {rng}" + ("" if ok else "   -> p is OUTSIDE this range, accuracy may drop"))
    print("  " + ftxt)
    print(THIN)
    show_path(xs, t, h, key, terms, show_table)
    print(THIN)
    for lab, c, d, _cells in terms:
        print(f"   {lab:<26} coef {fmt(c):>12}  x  {fmt(d):>12}  =  {fmt(c * d):>13}")
    print("   " + "-" * 66)
    print(f"   {'TOTAL':<26}{'':>34}  =  {fmt(total, 10):>13}")
    print(f"  RESULT ({name}): f({fmt(x)}) = {fmt(total, 10)}")


def show_central(xs, t, h, x, keys, cx, detail=True):
    n = len(xs)
    print("\n" + LINE)
    print(f"  QUERY  x = {fmt(x)}")
    print(LINE)
    lo = clamp(floor((x - xs[0]) / h), 0, n - 2)
    print(f"  Table range : [{fmt(xs[0])}, {fmt(xs[-1])}],  h = {fmt(h)}")
    print(f"  Query lies between {fmt(xs[lo])} and {fmt(xs[lo + 1])}")
    results = []
    for k in keys:
        fn, ftxt, rng, okf = CENTRAL[k]
        name, i, p, terms = fn(xs, t, h, x)
        total = sum(c * d for _l, c, d, _c in terms)
        res = (name, i, p, terms, total, ftxt, rng, okf(p))
        results.append(res)
        if detail:
            show_central_detail(xs, t, h, x, k, res, len(keys) == 1)
    if len(keys) == 1:
        print(LINE)
        return
    ref = poly_eval(cx, x)
    print("\n" + LINE)
    print(f"  COMPARISON at x = {fmt(x)}")
    print(LINE)
    print(f"  {'Formula':<16}{'origin':>8}{'p':>10}{'terms':>7}{'max δ':>7}{'result':>16}{'|err vs poly|':>15}  p-range")
    print("  " + "-" * 76)
    for name, i, p, terms, total, _, _, ok in results:
        print(f"  {name:<16}{fmt(xs[i]):>8}{fmt(p, 5):>10}{len(terms):>7}{max_order(terms):>7}"
              f"{fmt(total, 10):>16}{abs(float(total - ref)):>15.3e}  {'ok' if ok else 'outside'}")
    print("  " + "-" * 76)
    print(f"  {'Full polynomial':<16}{'':>8}{'':>10}{n:>7}{n - 1:>7}{fmt(ref, 10):>16}   (degree {len(cx) - 1}, all {n} points)")
    vals = [float(r[4]) for r in results]
    print(f"  Spread (max - min) over the {len(results)} formulae = {max(vals) - min(vals):.3e}")
    best = min(results, key=lambda r: (0 if r[7] else 1, abs(float(r[4] - ref))))
    print(f"  Best-suited formula here (p in its range, closest to full polynomial): {best[0]}")
    print(LINE)


# ================================================================== polynomial
def mul_linear(poly, root):
    out = [Fraction(0)] * (len(poly) + 1)
    for i, c in enumerate(poly):
        out[i + 1] += c
        out[i] -= c * root
    return out


def trim(c):
    while len(c) > 1 and abs(float(c[-1])) < 1e-12:
        c.pop()
    return c


def poly_in_p(t):
    n = len(t)
    tot = [Fraction(0)] * n
    basis = [Fraction(1)]
    for k in range(n):
        c = t[0][k] / factorial(k)
        for j, b in enumerate(basis):
            tot[j] += c * b
        basis = mul_linear(basis, k)
    return trim(tot)


def poly_in_x(xs, t, h):
    n = len(xs)
    tot = [Fraction(0)] * n
    basis = [Fraction(1)]
    for k in range(n):
        c = t[0][k] / (factorial(k) * h ** k)
        for j, b in enumerate(basis):
            tot[j] += c * b
        basis = mul_linear(basis, xs[0] + k * h)
    return trim(tot)


def poly_eval(c, x):
    r = Fraction(0)
    for a in reversed(c):
        r = r * x + a
    return r


def poly_str(c, var):
    parts = []
    for pw in range(len(c) - 1, -1, -1):
        a = c[pw]
        if abs(float(a)) < 1e-12:
            continue
        m = abs(a)
        s = "" if (m == 1 and pw > 0) else fmt(m)
        s += var if pw == 1 else (f"{var}^{pw}" if pw > 1 else "")
        parts.append(("-" if a < 0 else "+", s))
    if not parts:
        return "0"
    return ("-" if parts[0][0] == "-" else "") + parts[0][1] + \
           "".join(f" {sg} {s}" for sg, s in parts[1:])


def show_polynomial(xs, t, h):
    n = len(xs)
    cp, cx = poly_in_p(t), poly_in_x(xs, t, h)
    print("\n" + LINE)
    print("  MAXIMUM DEGREE POLYNOMIAL THROUGH THE DATA")
    print(LINE)
    print(f"  Maximum degree possible  : {n - 1}   (n = {n} points)")
    print(f"  Actual degree of the data: {len(cx) - 1}")
    print(f"\n  in terms of p  [ p = (x - {fmt(xs[0])}) / {fmt(h)} ] :")
    print(f"     f(p) = {poly_str(cp, 'p')}")
    print(f"\n  in terms of x :")
    print(f"     P(x) = {poly_str(cx, 'x')}")
    print(LINE)


# ================================================================== menus
def ask_queries(xs, action):
    """repeatedly ask for x; blank returns to the menu.  Only x INSIDE the table is accepted."""
    print(f"  (x must lie inside the table range [{fmt(xs[0])}, {fmt(xs[-1])}] - extrapolation is not supported)")
    while True:
        s = input("\nQuery x (blank = back to menu) : ").strip()
        if not s:
            return
        try:
            q = num(s)
        except (ValueError, ZeroDivisionError):
            print("  ! Invalid number.")
            continue
        if q < xs[0] or q > xs[-1]:
            print(f"  ! x = {fmt(q)} is outside [{fmt(xs[0])}, {fmt(xs[-1])}]  (that would be extrapolation) - "
                  "please enter a value inside the table.")
            continue
        action(q)


def central_menu(xs, t, h):
    cx = poly_in_x(xs, t, h)
    print_central_table(xs, t)
    while True:
        print("\n" + LINE)
        print("  CENTRAL DIFFERENCE FORMULAE")
        print(LINE)
        print("   1. Stirling's formula")
        print("   2. Bessel's formula")
        print("   3. Everett's formula")
        print("   4. Gauss Forward formula")
        print("   5. Gauss Backward formula")
        print("   6. Compare ALL central formulae")
        print("   0. Back to main menu")
        ch = input("  Choose formula : ").strip()
        if ch == "0":
            return
        if ch in CENTRAL:
            keys, detail = [ch], True
        elif ch == "6":
            keys = ["1", "2", "3", "4", "5"]
            detail = input("  Show term-by-term working for each formula? (y/n) [y] : ").strip().lower() != "n"
        else:
            print("  ! Invalid choice.")
            continue
        ask_queries(xs, lambda q: show_central(xs, t, h, q, keys, cx, detail))


def main():
    print(LINE)
    print("   NUMERICAL INTERPOLATION  -  NEWTON & CENTRAL DIFFERENCE FORMULAE")
    print(LINE)
    while True:
        try:
            n = int(input("Enter number of data points : "))
            if n >= 2:
                break
        except ValueError:
            pass
        print("  ! Enter an integer >= 2.")

    xs = read_list("Enter x value", n)
    if not equally_spaced(xs):
        print("\nError: x values must be distinct and equally spaced.")
        return

    fs = read_list("Enter f(x) value", n)
    pts = sorted(zip(xs, fs))
    xs, fs = [a for a, _ in pts], [b for _, b in pts]
    h = xs[1] - xs[0]
    if not equally_spaced(xs):
        print("\nError: x values must be distinct and equally spaced.")
        return

    t = build_table(fs)
    print_table(xs, t)

    while True:
        print("\n" + LINE)
        print("  MAIN MENU")
        print(LINE)
        print("   1. Newton FORWARD difference interpolation")
        print("   2. Newton BACKWARD difference interpolation")
        print("   3. CENTRAL difference interpolation  (Stirling / Bessel / Everett / Gauss)")
        print("   4. Show maximum-degree polynomial through the data")
        print("   5. Show difference table again")
        print("   0. Exit")
        ch = input("  Your choice : ").strip()
        if ch == "1":
            ask_queries(xs, lambda q: show_query(xs, t, h, q, "FORWARD"))
        elif ch == "2":
            ask_queries(xs, lambda q: show_query(xs, t, h, q, "BACKWARD"))
        elif ch == "3":
            central_menu(xs, t, h)
        elif ch == "4":
            show_polynomial(xs, t, h)
        elif ch == "5":
            print_table(xs, t)
            print_central_table(xs, t)
        elif ch == "0":
            print("\nGoodbye!")
            break
        else:
            print("  ! Invalid choice, enter 0-5.")


if __name__ == "__main__":
    main()