#running hospitalization calculation in parallel
nohup R CMD BATCH --no-save '--args minFacID=1 maxFacID=50' SyntheticSimulated.r synLog1.Rout &
nohup R CMD BATCH --no-save '--args minFacID=51 maxFacID=100' SyntheticSimulated.r synLog2.Rout &
nohup R CMD BATCH --no-save '--args minFacID=101 maxFacID=150' SyntheticSimulated.r synLog3.Rout &
nohup R CMD BATCH --no-save '--args minFacID=151 maxFacID=200' SyntheticSimulated.r synLog4.Rout &
nohup R CMD BATCH --no-save '--args minFacID=201 maxFacID=250' SyntheticSimulated.r synLog5.Rout &
nohup R CMD BATCH --no-save '--args minFacID=251 maxFacID=300' SyntheticSimulated.r synLog6.Rout &
nohup R CMD BATCH --no-save '--args minFacID=301 maxFacID=350' SyntheticSimulated.r synLog7.Rout &
nohup R CMD BATCH --no-save '--args minFacID=351 maxFacID=400' SyntheticSimulated.r synLog8.Rout &
nohup R CMD BATCH --no-save '--args minFacID=401 maxFacID=450' SyntheticSimulated.r synLog9.Rout &
#running mortality calculation in parallel
nohup R CMD BATCH --no-save '--args minFacID=1 maxFacID=50' SyntheticSimulatedMort.r synLogD1.Rout &
nohup R CMD BATCH --no-save '--args minFacID=51 maxFacID=100' SyntheticSimulatedMort.r synLogD2.Rout &
nohup R CMD BATCH --no-save '--args minFacID=101 maxFacID=150' SyntheticSimulatedMort.r synLogD3.Rout &
nohup R CMD BATCH --no-save '--args minFacID=151 maxFacID=200' SyntheticSimulatedMort.r synLogD4.Rout &
nohup R CMD BATCH --no-save '--args minFacID=201 maxFacID=250' SyntheticSimulatedMort.r synLogD5.Rout &
nohup R CMD BATCH --no-save '--args minFacID=251 maxFacID=300' SyntheticSimulatedMort.r synLogD6.Rout &
nohup R CMD BATCH --no-save '--args minFacID=301 maxFacID=350' SyntheticSimulatedMort.r synLogD7.Rout &
nohup R CMD BATCH --no-save '--args minFacID=351 maxFacID=400' SyntheticSimulatedMort.r synLogD8.Rout &
nohup R CMD BATCH --no-save '--args minFacID=401 maxFacID=450' SyntheticSimulatedMort.r synLogD9.Rout &
