set more off
use ./3_Temp/RStataDataInmisspecify_SEIR647
capture noisily {
/*RSTATA: cut me here*/
set matsize 1000
      glm Rt i.unit i.week 1.trt_post, family(poisson) link(log)
      boottest 1.trt_post, cluster(unit) reps(10000) quietly
      gen p = r(p) in 1
      keep p
      keep if _n==1
/*RSTATA: cut me here*/
} /* end capture noisily */
saveold ./3_Temp/RStataDataOutmisspecify_SEIR647, version(12)
exit, clear STATA
