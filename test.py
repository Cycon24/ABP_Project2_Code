from freeVortexCompressor import compressorStageData, freeVortexCompressorMeanLine

obj = compressorStageData({"M_1" : 1}, "test")
data = obj.toDataFrame(1)
print(data)