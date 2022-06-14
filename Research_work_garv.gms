$title Mathematical Model for Cattle Feed Mix Problem With Chance Constraint
$onText
G.Singhal1 , M.Mann2
IIIT Sonepat
Reasearch Work
$offText
Set
   f 'feeds'     / Grass, Wheat-Straw, Mustard-Cake, Wheat-bran, Binola-Khal, Chickpea-peel, Porridge /
   n 'nutrients' / protein, fiber /;
Parameter
   price(f) 'feed prices (feed per kg)' /  Wheat-Straw    170.00
                                           Mustard-Cake    30.00
                                           Wheat-bran       7.00
                                           Binola-khal     39.00
                                           Chickpea-peel   90.00
                                           Porridge        18.00 /
   req(n)   'requirements         (pct)' / protein   10.7
                                           fiber    22.0/;

Table char(*,n,f) 'feed characteristics (pct)'
                      Wheat-Straw  Mustard-Cake  wheat-bran Binola-khal Chickpea-peel porridge
   mean.protein        6.0          36.0           12.7        22.5        8.86       11.22
   mean.fiber         46.2          12.0           16.9        27.4       22.4        15.00
   variance.protein     .28           .19          20.5          .62       1.2         2.3
   variance.fiber      2.34           1.19         10.5          2.2       2.2        11.5;

Variable
   cost 'total cost per kg'
   x(f) 'feed mix     (pct)';

Positive Variable x;

Equation
   cdef    'cost definition'
   mc      'mix constraint'
   nbal(n) 'nutrient balance'
   cc(n)   'chance constraint';

cdef..    cost =e= sum(f, price(f)*x(f));

mc..      sum(f, x(f)) =e= 1;

nbal(n).. sum(f, char("mean",n,f)*x(f)) =g= req(n);

cc(n)..   sum(f, char("mean",n,f)*x(f)) - 1.645*sqrt(sum(f, char("variance",n,f)*sqr(x(f)))) =g= req(n);
Model
   det    'deterministic model' / cdef, mc, nbal    /
   chance 'chance model'        / cdef, mc,      cc /;

Parameter mix 'mixing report';

option lp = cplex;
solve det    minimizing cost using lp;
mix(f,'det   ') = x.l(f);

option nlp = conopt4;
solve chance minimizing cost using nlp;
mix(f,'chance') = x.l(f);

display mix;