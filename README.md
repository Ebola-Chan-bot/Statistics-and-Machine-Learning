埃博拉酱的统计和机器学习工具箱工具箱，提供一系列 MATLAB Statistics and Machine Learning Toolbox 所欠缺，但却常用的增强功能。

个别功能最高要求 MATLAB R2025a 版本才能正常工作，但低版本一般可以运行大部分功能。

本工具箱中所有函数均在`StatisticsAndMachineLearning`命名空间下，使用前需import。

```MATLAB
%对多维数组沿指定维度执行PCA分析
function [Coeff,Score,Explained] = DimensionalPca(Array,Dimensions,NumComponents)

%内置anovan的升级版，将变量表作为分组输入，额外支持多重比较
function varargout = TabularAnovaN(varargin)
```