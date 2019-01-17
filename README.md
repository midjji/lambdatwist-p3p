# lambdatwist-p3p 
A simple, easy to understand p3p solver with state of the art performance.
http://openaccess.thecvf.com/content_ECCV_2018/html/Mikael_Persson_Lambda_Twist_An_ECCV_2018_paper.html


first cleaned draft,


todo
verfy shared seed and generator clang and gcc between various versikns 
slowdown on ubuntu 18 in some config?  unverified, chck. 

copy docs from old benchmark

rewrite without dependency in cvl for opencv and/or opengv. 
add arxiv version with the fixed eqn 3. 

add poster? fix eqn 3.

# Duplicate solutions
Its a bit tricky to distinguish duplicate solutions from nearby solutions,a threshold 10^-5 is used.

# Benchmark results!
Im working on a website


Full tests - computer i8
-----------------------------------------------------------------------------
|  Table                                    | Lambda   | Ke       | Kneip    |
|  ground truth found                       | 9999980  | 9998929  | 9995835  |
|  any solution found                       | 9999996  | 9999771  | 9999333  |
|  ground truth found ratio                 | 0.999998 | 0.999893 | 0.999583 |
|  no solution found                        | 4        | 229      | 667      |
|  valid according to solver                | 16883653 | 17076564 | 16882854 |
|  duplicates                               | 1        | 191453   | 370      |
|  unique correct solutions                 | 16883646 | 16885111 | 16882484 |
|  incorrect solutuons output by the solver | 6        | 0        | 0        |
-----------------------------------------------------------------------------

---------------------------------------------------------------------------
| Timer        | Total   | Mean    | Median  | Min     | Max     | Samples |
| Lambda Total | 2102ms  | 2102ms  | 2102ms  | 2102ms  | 2102ms  | 1       |
| ke Total     | 3634ms  | 3634ms  | 3634ms  | 3634ms  | 3634ms  | 1       |
| Kneip Total  | 9965ms  | 9965ms  | 9965ms  | 9965ms  | 9965ms  | 1       |
---------------------------------------------------------------------------


Minimal tests, the other algorithms are faster but get bad results:
On computer - computer i8
  -----------------------------------------------------------------------------
 |  Table                                    | Lambda   | Ke       | Kneip    |
 |  ground truth found                       | 9999980  | 9998969  | 9995835  |
 |  any solution found                       | 9999996  | 9999771  | 9999333  |
 |  ground truth found ratio                 | 0.999998 | 0.999897 | 0.999583 |
 |  no solution found                        | 4        | 229      | 667      |
 |  valid according to solver                | 16883653 | 26706936 | 24240084 |
 |  duplicates                               | 1        | 191453   | 370      |
 |  unique correct solutions                 | 16883646 | 16885111 | 16882484 |
 |  incorrect solutuons output by the solver | 6        | 9630372  | 7357230  |
  -----------------------------------------------------------------------------

---------------------------------------------------------------------------
| Timer        | Total   | Mean    | Median  | Min     | Max     | Samples |
| Lambda Total | 2098ms  | 2098ms  | 2098ms  | 2098ms  | 2098ms  | 1       |
| ke Total     | 2586ms  | 2586ms  | 2586ms  | 2586ms  | 2586ms  | 1       |
| Kneip Total  | 8489ms  | 8489ms  | 8489ms  | 8489ms  | 8489ms  | 1       |










