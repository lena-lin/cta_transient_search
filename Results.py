import numpy as np

alarm_counts = np.array([1,0,0,0,1,0,1,2,1,0,2,1,1,1,1,2,0,1,0,2,0,0,0,
                         2,1,2,0,3,1,0,0,2,1,0,1,1,1,2,1,2,0,1,2,0,1,3,
                         3,1,2,1,1,0,0,1,0,2,0,1,1,1,0,2,1,0,0,1,1,1,0,
                         0,0,0,1,2,0,1,0,2,0,0,0,0,2,2,1,0,0,1,0,1,4,0,
			 1,1,2,2,1,1,0,0,])

Average = np.mean(alarm_counts)
Uncertainties = np.std(alarm_counts)
print('Simulated ',len(alarm_counts) ,'times 8 h of observation with each 48 time series')
print('Averaged 8 h of observation with ', Average, ' pm ', Uncertainties, 'False alarms')
print('With 1000h of observation per year, calculate 125 8h slices, resulting in ' , Average*125, 'pm', Uncertainties*125 ,' Flase alarms')
