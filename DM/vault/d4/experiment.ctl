method ml
set szr 0
set sv 0
set st0 0
set p 0
depends v target1
depends v target2
depends v lag
depends a target1
depends a target2
depends a lag
depends t0 target1
depends t0 target2
depends t0 lag
format RESPONSE TIME target1 target2 lag
load *.txt
log blinkFinal.txt