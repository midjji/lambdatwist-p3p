# lambdatwist-p3p 
A simple, easy to understand p3p solver with state of the art performance.
http://openaccess.thecvf.com/content_ECCV_2018/html/Mikael_Persson_Lambda_Twist_An_ECCV_2018_paper.html


This is the benchmark code. 

A corresponding pnp solver is found in:
https://github.com/midjji/pnp

For a wrapper around the p3p specifically look to:
https://github.com/midjji/pnp/blob/master/p4p.h

 |-------------------------------------------|----------|----------|----------|
 |  Method                                   | Lambda   | Ke       | Kneip    |
 |  ground truth found                       | 9999991  | 9997237  | 9990804  |
 |  any solution found                       | 9999996  | 9999663  | 9999228  |
 |  ground truth found ratio                 | 0.999999 | 0.999724 | 0.99908  |
 |  no solution found                        | 4        | 337      | 772      |
 |  valid according to solver                | 16827636 | 17811966 | 24152163 |
 |  duplicates                               | 0        | 163697   | 379      |
 |  unique correct solutions                 | 16827635 | 16852892 | 16826255 |
 |  incorrect solutuons output by the solver | 0        | 795377   | 7325529  |
 |-------------------------------------------|----------|----------|----------|
 
|--------------|---------|---------|---------|---------|---------|----------|
| Timer        | Total   | Mean    | Median  | Min     | Max     | Samples | 
| Lambda Total | 12356ms | 2471ms  | 2469ms  | 2467ms  | 2478ms  | 5       | 
| ke Total     | 17241ms | 3448ms  | 3449ms  | 3436ms  | 3465ms  | 5       | 
| Kneip Total  | 29s     | 5887ms  | 5870ms  | 5863ms  | 5965ms  | 5       | 
|--------------|---------|---------|---------|---------|---------|----------|
