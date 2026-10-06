# Viva Notes: Interpolation with Newton and Central Differences

Covers Newton Forward, Newton Backward, Gauss Forward, Gauss Backward, Stirling, Bessel and Everett.

---

## 1. Basics

**Interpolation** means finding f(x) for an x that lies *inside* the table. Estimating a value outside the table is *extrapolation*. My program only does interpolation and rejects x outside [x_first, x_last].

All seven formulae need **equally spaced** x (step h). All use the same variable:

> **p = (x − x₀) / h**, so x = x₀ + p·h

p is the number of steps the query is from the origin x₀. It is a distance, not a table index.

### Difference operators

| Operator | Definition | Where it points |
|---|---|---|
| Forward Δ | Δfᵣ = fᵣ₊₁ − fᵣ | down-right |
| Backward ∇ | ∇fᵣ = fᵣ − fᵣ₋₁ | up-right |
| Central δ | δfᵣ = f(r+½) − f(r−½) | symmetric about r |
| Shift E | E fᵣ = fᵣ₊₁ | E = 1 + Δ = (1 − ∇)⁻¹ |

All three difference tables hold the same numbers. Only the names (the subscripts) change:

> Δᵏf_r = ∇ᵏf_(r+k) = δᵏf at the half-way position r + k/2

So δ²f₀ = Δ²f₋₁ and δ³f_(1/2) = Δ³f₋₁.

**Rule of thumb:** even-order central differences (δ², δ⁴) sit on the same row as an x-value. Odd-order ones (δ, δ³, δ⁵) sit between two rows, at half-integer positions.

---

## 2. The formulae, paths and ranges

| Formula | Origin x₀ | Range of p | Differences used | Path through the table |
|---|---|---|---|---|
| **Newton Forward** | first point of the part you use (point just before x) | 0 ≤ p < 1 | Δ at x₀ | straight **down** the diagonal |
| **Newton Backward** | last point of the part you use | 0 ≤ p < 1 (or −1 < p ≤ 0 when taken from the point after x) | ∇ at x₀ | straight **up** the diagonal |
| **Gauss Forward** | point just before x | 0 ≤ p < 1 | δ at x₀, alternating | zig-zag, goes **down** first |
| **Gauss Backward** | point just after x | −1 < p ≤ 0 | δ at x₀, alternating | zig-zag, goes **up** first |
| **Stirling** | **nearest** point to x | −¼ < p < ¼ | **mean** of the two odd δ, even δ at x₀ | zig-zag through the **middle row** |
| **Bessel** | point just before x | p near ½ (¼ to ¾) | odd δ at x₀+½, **mean** of two even δ | zig-zag **between two rows** x₀ and x₁ |
| **Everett** | point just before x | 0 ≤ p < 1 | **even differences only** on rows x₀ and x₁ | **two horizontal lines** |

### Formulae

```
Newton Fwd : f = f0 + pΔf0 + p(p-1)/2! Δ²f0 + p(p-1)(p-2)/3! Δ³f0 + ...
Newton Bwd : f = f0 + p∇f0 + p(p+1)/2! ∇²f0 + p(p+1)(p+2)/3! ∇³f0 + ...

Gauss Fwd  : f = f0 + pδf(1/2) + p(p-1)/2! δ²f0 + (p+1)p(p-1)/3! δ³f(1/2)
                 + (p+1)p(p-1)(p-2)/4! δ⁴f0 + ...
Gauss Bwd  : f = f0 + pδf(-1/2) + (p+1)p/2! δ²f0 + (p+1)p(p-1)/3! δ³f(-1/2)
                 + (p+2)(p+1)p(p-1)/4! δ⁴f0 + ...

Stirling   : f = f0 + p/2 [δf(-1/2)+δf(1/2)] + p²/2! δ²f0
                 + p(p²-1)/3! · ½[δ³f(-1/2)+δ³f(1/2)] + p²(p²-1)/4! δ⁴f0 + ...

Bessel     : f = f0 + pδf(1/2) + p(p-1)/2! · ½[δ²f0+δ²f1]
                 + p(p-1)(p-½)/3! δ³f(1/2) + (p+1)p(p-1)(p-2)/4! · ½[δ⁴f0+δ⁴f1] + ...

Everett    : f = q f0 + q(q²-1)/3! δ²f0 + q(q²-1)(q²-4)/5! δ⁴f0 + ...
               + p f1 + p(p²-1)/3! δ²f1 + p(p²-1)(p²-4)/5! δ⁴f1 + ...       (q = 1 − p)
```

---

## 3. How the formulae are related

1. **Stirling = mean of Gauss Forward and Gauss Backward**, both taken from the same origin x₀. This is why Stirling uses *means* of two odd differences.
2. **Bessel = mean of Gauss Forward (origin x₀) and Gauss Backward (origin x₁)**. It is symmetric about the middle of the interval [x₀, x₁], which is why it works best near p = ½.
3. **Everett** comes from Bessel/Gauss by replacing every odd difference using δ^(2m+1)f(½) = δ^(2m)f₁ − δ^(2m)f₀. After that only even differences remain. It is Lagrange-like: f = (q-part from row x₀) + (p-part from row x₁).
4. **Newton Forward and Backward** are the "one-sided" versions. They take all differences from one diagonal, not from the centre.

**Key point:** if all of them use the same set of data points and go to the full order, they give the **same polynomial and the same answer**. In my program, the full-degree polynomial through all points is printed as a reference. The formulae only differ in:
- how fast the terms shrink (accuracy when truncated),
- which table entries they need, and
- ease of calculation.

---

## 4. Why p-ranges and accuracy differ

- **Stirling, p small (|p| < ¼).** The term coefficients contain factors p², p²−1, p²−4, so they are very small when p is near 0. It is most accurate close to a tabular point.
- **Bessel, p near ½.** The 3rd-difference coefficient has the factor **(p − ½)**. It is **exactly 0 at p = ½**, and all odd-order terms of the form (p − ½)·(...) vanish at the mid-point. So the odd differences contribute nothing and the series converges fastest at the centre of the interval. Bessel is the most popular formula.
- **Everett, any p in [0,1).** It needs only even differences, so you do not compute odd columns. The coefficients depend only on p (for the p-part) and q = 1 − p (for the q-part), so they are easy to tabulate. It is simple and fast.
- **Gauss Forward / Backward.** Valid right at x₀ with p in [0,1) or (−1,0]. They are the "parents" of Stirling and Bessel. They are less accurate than the means because the error terms of the two Gauss formulae partly cancel when averaged.
- **Newton Forward / Backward.** Best at the very **beginning** or **end** of the table, where no central differences exist (no points on both sides). Near the centre they use fewer, less balanced differences than the central formulae.

### Which formula do I pick?

| Where x lies | Best choice |
|---|---|
| In the first interval (start of table) | Newton Forward |
| In the last interval (end of table) | Newton Backward |
| Middle of table, very close to a tabular point (|p| < ¼) | **Stirling** |
| Middle of table, about halfway between points (p ≈ ½) | **Bessel** |
| Middle of table, only even differences wanted | **Everett** |
| Middle, but one-sided path wanted | Gauss Forward (p ≥ 0) / Gauss Backward (p ≤ 0) |

---

## 5. Likely viva questions and answers

**Q1. What is the difference between interpolation and extrapolation?**
Interpolation estimates f(x) inside the table range. Extrapolation estimates it outside the range and is less reliable. My program only does interpolation, so an x outside the table is rejected.

**Q2. Why must the data be equally spaced?**
All these formulae are built from finite differences and the shift operator E = 1 + Δ, which need a constant step h so that x = x₀ + ph. For unequal spacing use Lagrange, Newton divided differences or Aitken.

**Q3. What is p?**
p = (x − x₀)/h, the number of steps from the origin to the query.

**Q4. Forward and backward differences: what is the difference?**
Δfᵣ = fᵣ₊₁ − fᵣ and ∇fᵣ = fᵣ − fᵣ₋₁. They are the same numbers in the same table, just named from different positions: Δᵏf_r = ∇ᵏf_(r+k).

**Q5. Why does Newton Backward have (p+1), (p+2)… and Forward have (p−1), (p−2)…?**
Forward uses E = 1 + Δ, so Eᵖ = (1+Δ)ᵖ gives binomial coefficients p(p−1)…. Backward uses E = (1−∇)⁻¹, so Eᵖ = (1−∇)⁻ᵖ gives p(p+1)….

**Q6. What is a central difference?**
δfᵣ = f(r+½) − f(r−½). δ²fᵣ = f(r+1) − 2fᵣ + f(r−1) = Δ²f(r−1). Even orders fall on the same row as x₀. Odd orders fall between rows.

**Q7. Why does Stirling use means?**
The odd central differences δᵏf(−½) and δᵏf(½) lie on either side of x₀. Averaging them puts the term on the middle row, and it makes Stirling the mean of Gauss Forward and Gauss Backward.

**Q8. Why is Stirling only good for small p?**
Because it is centred on x₀. The coefficients p, p², p(p²−1)… are smallest when p is near 0, so the series converges fastest there and the formula is recommended for −¼ < p < ¼.

**Q9. Why is Bessel best near p = ½?**
It is centred between x₀ and x₁. Its odd-difference coefficients contain (p − ½) and vanish at p = ½. It is also the mean of Gauss Forward from x₀ and Gauss Backward from x₁.

**Q10. What is special about Everett?**
It uses only even differences (δ², δ⁴, …) along rows x₀ and x₁, with q = 1 − p. It is easy to program and fast, and it is used for double interpolation in tables.

**Q11. Do different formulae give different answers?**
If they use all the same data and go to full order they give the same value (the same polynomial). Differences arise from **truncation**: each formula stops after a different number of terms or uses different differences, so their errors are different. Near the table edges central formulae run out of differences and lose accuracy.

**Q12. How do you read a central difference from the forward table?**
δ^(2m)f_i = Δ^(2m)f_(i−m) and δ^(2m+1)f_(i+½) = Δ^(2m+1)f_(i−m).

**Q13. How many points are needed, and what degree can you reach?**
With n points the highest difference is order n − 1, so the highest polynomial degree is n − 1. If the k-th differences are constant (and the next are zero), the data is exactly a degree-k polynomial.

**Q14. Why does the program use fractions?**
Exact arithmetic means no rounding error, so any small difference between formulae is purely truncation and not floating-point noise.

**Q15. Why would the Gauss Backward δ⁴ coefficient be (p+2)(p+1)p(p−1)/4!?**
It follows the pattern C(p + ⌊k/2⌋, k). The manual's version for this term looks like the forward one and is likely a typo.

**Q16. What if p is outside the recommended range?**
The formulae are still algebraically valid, but more terms are needed for the same accuracy. The program flags this ("p-range: outside") in the comparison table.

---

## 6. One-line summary of each formula

- **Newton Forward:** one diagonal downward from x₀. Use it at the **start** of the table.
- **Newton Backward:** one diagonal upward to x₀. Use it at the **end** of the table.
- **Gauss Forward:** zig-zag with the first step down. A central formula for p in [0,1).
- **Gauss Backward:** zig-zag with the first step up. A central formula for p in (−1,0].
- **Stirling:** the average of Gauss F and B. Use it when x is **close to a tabular point** (|p| < ¼).
- **Bessel:** works between two rows. Use it when x is **halfway between points** (p ≈ ½).
- **Everett:** **even differences only** from the two surrounding rows. Simple and fast.

---

## 7. How the program shows the path

Each formula prints a **PATH FOLLOWED** section, as in the lab manual. It lists the table entries used, in order, with their positions and values. The difference table is then reprinted with those entries marked in `[ ]`, so you can see the zig-zag, diagonal or horizontal lines directly on the table.