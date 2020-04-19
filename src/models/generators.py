import memory_profiler
import time
import pandas as pd

csv_gen = pd.read_csv("data/model_input.csv", sep=",", header=0, parse_dates=["When"])

# csv_gen = (row for row in open("data/model_input.csv"))

# df = pd.read_csv("data/model_input.csv", sep=",", header=0, parse_dates=["When"])


m1 = memory_profiler.memory_usage()
t1 = time.time()

row_count = 0

for row in csv_gen.itertuples():
    print(row.When)
print(f"Row count is {row_count}")

t2 = time.time()
m2 = memory_profiler.memory_usage()
time_diff = t2 - t1
mem_diff = m2[0] - m1[0]
print(f"It took {time_diff} Secs and {mem_diff} Mb to execute this method")