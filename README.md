# lambdatwist-p3p 
A simple, easy to understand p3p solver with state of the art performance.
http://openaccess.thecvf.com/content_ECCV_2018/html/Mikael_Persson_Lambda_Twist_An_ECCV_2018_paper.html

As of 2022 the solver has been patched to address a missing special case which is described in my thesis: 
https://www.diva-portal.org/smash/record.jsf?pid=diva2%3A1635583&dswid=3393

Cite the latter if the special case is relevant, otherwise cite the ECCV publication. Preferably cite both of course :P. 



This is the benchmark code. 

A corresponding pnp solver is found in:
https://github.com/midjji/pnp

For a wrapper around the p3p specifically look to:
https://github.com/midjji/pnp/blob/master/p4p.h

 
 |  Method                                   | Lambda   | Ke       | Kneip    |
 |-------------------------------------------|----------|----------|----------|
 |  ground truth found                       | 9999991  | 9997237  | 9990804  |
 |  any solution found                       | 9999996  | 9999663  | 9999228  |
 |  ground truth found ratio                 | 0.999999 | 0.999724 | 0.99908  |
 |  no solution found                        | 4        | 337      | 772      |
 |  valid according to solver                | 16827636 | 17811966 | 24152163 |
 |  incorrect solutions output by the solver | 0        | 795377   | 7325529  |

Lambdatwist provides all unique correct solutions, solutions output by Ke, or Kneip in addition to these are either not unique, or not correct. 


| Timer   | Total   | Mean    | Median  | Max     | Samples |
|---------|---------|---------|---------|---------|---------|
| Lambda  | 12s     | 2471ns  | 2469ns  | 2478ns  | 5000000 |
| Ke      | 17s     | 3448ns  | 3449ns  | 3465ns  | 5000000 | 
| Kneip   | 29s     | 5887ns  | 5870ns  | 5965ns  | 5000000 |

