import nidaqmx
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# Read from DAQ Device
def readdaq():
    task = nidaqmx.Task()
    task.ai_channels.add_ai_voltage_chan("Dev1/ai9")
    task.start()
    value = task.read()
    task.stop()
    task.close()
    return value
# Write Data Function
def writefiledata(t, x):
# Open File
    file = open("tempdata.txt", "a")
    # Write Data
    time = str(t)
    value = str(round(x, 2))
    file.write(time + "\t" + value)
    file.write("\n")
    # Close File
    file.close()
# Initialize Logging
Ts = 1 # Sampling Time [seconds]
N = 10
k = 1
x_len = N # Number of points to display
Tmin = -2; Tmax = 2
y_range = [Tmin, Tmax] # Range of possible Y values to display
data = []
# Create figure for plotting
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
xs = list(range(0, N))
ys = [0] * x_len
ax.set_ylim(y_range)

# Create a blank line. We will update the line in animate
line, = ax.plot(xs, ys, '.')
# Configure Plot
plt.title('Voltage')
plt.xlabel('t [s]')
plt.ylabel('Voltage [volt]')
plt.grid()
#Logging Temperature Data from DAQ Device
def logging(i, ys):
    value = readdaq()
    print("V =", round(value,1), "[volts]")
    data.append(value)
    time.sleep(Ts)
    global k
    k = k + 1
    # writefiledata(k*Ts, value)
    # Add y to list
    ys.append(value)
    # Limit y list to set number of items
    ys = ys[-x_len:]
    # Update line with new Y values
    line.set_ydata(ys)
    return line,
ani = animation.FuncAnimation(fig,logging,fargs=(ys,),interval=100,blit=True)
# plt.show()