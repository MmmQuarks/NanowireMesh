import time

class Timer:
    def __init__(self):
        self.times = {}

    # starts timer for process.
    def start(self, process):
        times = self.times
        if process in times.keys():
            if times[process]['start'] != 0:
                print(''.join(['Error starting timer for process: ', process, '. Previous timer not stopped.']))
            else:
                times[process]['start'] = time.time()
        else:
            times[process] = {'total' : 0, 'start' : time.time()}
        self.times = times

    def stop(self, process):
        times = self.times
        if process in times.keys():
            times[process]['total'] += time.time() - times[process]['start']
            times[process]['start'] = 0
        else:
            print(''.join(['Error stopping timer for process: ', process, '. No timer created.']))
        self.times = times
