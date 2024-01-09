import matplotlib.pyplot as plt

filename = 'result_mse.txt'
num = []
with open(filename, 'r') as file:
    while True:
        line = file.readline()
        if not line:
            break
        line = line[-9:]
        num.append(float(line) / 10)

n = 16
t = range(-5, 21)
plt.plot(t[:n], num[:n])    
plt.xlabel('SNR(dB)')
plt.ylabel('RMSE')
plt.title('RMSE')

plt.savefig('data/RMSE.png')
plt.show()