#run examineRealization "/home/pwg/Realizations/nokrige.hdf5" 0 0 False TRUE
#run examineRealization "/home/pwg/Realizations/nokrige.hdf5" 1 0 False TRUE
#run examineRealization "/home/pwg/Realizations/nokrige.hdf5" 2 0 False TRUE
#run examineRealization "/home/pwg/Realizations/nokrige.hdf5" 3 0 False TRUE
#run examineRealization "/home/pwg/Realizations/nokrige.hdf5" 4 0 False TRUE
#run examineRealization "/home/pwg/Realizations/nokrige.hdf5" 5 0 False TRUE
#run examineRealization "/home/pwg/Realizations/nokrige.hdf5" 6 0 False TRUE
#run examineRealization "/home/pwg/Realizations/nokrige.hdf5" 7 0 False TRUE
#run examineRealization "/home/pwg/Realizations/nokrige.hdf5" 8 0 False TRUE

#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_6_7.hdf5" 0 270 True FALSE
#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_40_41.hdf5" 0 270 True FALSE
#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_235_236.hdf5" 0 270 True FALSE
#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_244_245.hdf5" 0 270 True FALSE
#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_251_252.hdf5" 0 270 True FALSE
#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_277_278.hdf5" 0 270 True FALSE
#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_416_417.hdf5" 0 270 True FALSE
#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_447_448.hdf5" 0 270 True FALSE
#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_496_497.hdf5" 0 270 True FALSE
#run examineRealization "/home/pwg/Realizations/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_503_504.hdf5" 0 270 True FALSE

#import time as time
#time.sleep(3600*12)

#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/krige.hdf5" 0 6 True TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/krige.hdf5" 1 6 True TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/krige.hdf5" 2 6 True TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/krige.hdf5" 3 6 True TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/krige.hdf5" 4 6 True TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/krige.hdf5" 5 6 True TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/krige.hdf5" 5 6 True TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/krige.hdf5" 6 6 True TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/krige.hdf5" 7 6 True TRUE

#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5" 0 6 False TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5" 1 6 False TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5" 2 6 False TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5" 3 6 False TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5" 4 6 False TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5" 5 6 False TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5" 5 6 False TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5" 6 6 False TRUE
#run examineRealization "/home/anand/research/mbg-world/mbgw-scripts/.nokrige.hdf5" 7 6 False TRUE

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5  2000 1.e8 5 2 1
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_1.hdf5" 0 1 False TRUE 1

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5  2000 1.e8 5 2 2
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_2.hdf5" 0 1 False TRUE 2

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5  2000 1.e8 5 2 3
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_3.hdf5" 0 1 False TRUE 3

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5  2000 1.e8 5 2 4
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_4.hdf5" 0 1 False TRUE 4

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5  2000 1.e8 5 2 5
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_5.hdf5" 0 1 False TRUE 5

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5  2000 1.e8 5 2 6
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_6.hdf5" 0 1 False TRUE 6

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5  2000 1.e8 5 2 7
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_7.hdf5" 0 1 False TRUE 7

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5  2000 1.e8 5 12 8
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_8.hdf5" 0 11 False TRUE 8

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5  2000 1.e8 5 2 9
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_9.hdf5" 0 1 False TRUE 9

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5 2000 1.e8 5 12 10 10000
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_10.hdf5" 0 11 False TRUE 10


cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5 2000 1.e8 5 12 11 10000
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_11.hdf5" 0 11 False TRUE 11


cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5 2000 1.e8 5 12 12 10000
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/nokrige-thick_12.hdf5" 0 11 False TRUE 12


cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5 2000 1.e8 5 13 10000
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/test_inAndout_krige_13.hdf5" 0 11 True TRUE 13

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5 2000 1.e8 5 14 10000
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/test_inAndout_krige_14.hdf5" 0 4 True TRUE 14

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5 2000 1.e8 5 15 15000
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/test_inAndout_krige_15.hdf5" 0 4 True TRUE 15

cd
cd /home/pwg/mbg-world/mbgw-scripts
run joint-sim-condor-zeus.py 0 1 1 AF ../datafiles/good-traces/QRYPFPR010708_Africa_Run_9.10.2008.hdf5 2000 1.e8 5 16 500
cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/test_inAndout_krige_16.hdf5" 0 4 True TRUE 16


#cd
#cd /home/pwg/mbg-world/mbgw-scripts
#run joint-sim-condor-zeus.py 0 1 1 AS ../datafiles/good-traces/QRYPFPR230708_Asia_Run_1.9.2008.hdf5  2000 1.e8 5 12 9 10000
#cd
#cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
#run examineRealization "/home/pwg/mbg-world/mbgw-scripts/AS_nokrige-thick_9.hdf5" 0 11 False TRUE 9


cd
cd /home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm
run examineRealization "/home/pwg/mbg-world/mbgw-scripts/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_0_1.hdf5" 0 4 True TRUE 15


cd
cd /home/pwg/mbg-world/mbgw-script
run examineRealization "/home/pwg/Realizations/postmapview/realizations_mem_100000000_QRYPFPR010708_Africa_Run_9.10.2008_iterations_0_1.hdf5" 0 4 True FALSE 15 0 11 True True

