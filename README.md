# Estimation Project #

### Step 1: Sensor Noise ###
This step involves adding noise to the sensor measurements to simulate a realistic environment.
To do this I took the graphs that were created from scenario 6 and put the data into two numpy arrays.
Then to find MeasuredStdDev_AccelXY and MeasuredStdDev_GPSPosXY I took the standard deviation of each array.

### Step 2: Attitude Estimation ###
TO estimate the Attitude of the drone the QuadEstimatorEKF uses the function UpdateFromIMU.

This function takes in the accerlation and gyro vectors and updates the roll estimation, pitch estimation, and gyro values.

To do this the function creats a rotation matrix based on the Euler angles to find the integration of the angles by multiplying the matrix by the gyro.

Then it uses these values to update the pitch and roll.


### Step 3: Prediction Step ###
This step involves setting up the prediction step of the estimator.
The function PredictState along with GetRbgPrime accomplishes this.

The GetRbgPrime function sets up the RbgPrime matrix which is the partial derivative of the Rbg rotation matrix.

The PredictState function uses the current state, time step, body fram acceleration, and gyro body rates to find the predicted state vector. To accomplish this the function takes respective current states and add the time step / integrator product.

### Step 4: Magnetometer Update ###
To havethe Magnetometer update correctly we use the UpdateFromMag function. This function takes in the yaw value of the Magnetometer and updates the EKF by calling the update function with the updated z, hPrime, magnitude of R, and the zFromX array.

It first takes the difference between the magnitude of the current yaw and the value of yaw in the state vector. 
Then it uses the difference to updates the zFromX vector and makes the hPrime value at 0, 6 = to 1.

### Step 5: Closed Loop + GPS Update ###
To get closed loop GPS update the function UpdateFromGPS takes in a position and velocity vector. zFromX is then updated to the values of the ekfState vector, and hPrime (0, 0), (1,1), (2,2), (3,3), (4,4), (5,5) are set to 1.
This function then takes these values and passes them into the Updatefunction; Which takes in z, hPrime, magnitude of R, and the zFromX array.


### Step 6: Adding Your Controller ###
For the final step I pulled in the code from my controller in the previous project and saw everything working as expected.