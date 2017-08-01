# trajOptimization
This project is for optimizing seperated trajactories data from crowd-sourcing sensor data. The data only contains consensus locations' reference in different trajectories.<br>
eg. <br>
>traj:7	lmGTID:10 1 2 9 7 4	lmIndex:11 44 65 99 124 157<br>
>traj:9	lmGTID:10 1 2 9 7 4	lmIndex:11 45 66 99 124 157<br>
>traj:10	lmGTID:10 1 2 9 7 4	lmIndex:11 45 66 99 124 157<br>
>traj:21	lmGTID:2 9 37 64 88 89	lmIndex:10 45 50 61 72 80<br>
>traj:25	lmGTID:2 9 37 64 88 89	lmIndex:10 44 49 60 71 80<br>
>(lmGTID refers to the same location' reference in global consesus location list.<br>
>lmIndex refers to these consensus location's refernce in each single trajectory.)<br>
     
# Dependency
g2o<br>
cholmod<br>
eigen<br>
  
# Compilation
Written in kdevelop (https://www.kdevelop.org/)<br>
Based on cmake.<br>
  
# Data
From Shenzhen University -》 Spatial Lab.<br>

