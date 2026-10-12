# 数据与复现范围

核心仓库保留 `formulations.csv`（样品浓度和训练/测试划分），不上传原始实验曲线和辅助处理数据。

## 使用已有成果

`artifacts/model_c/model_c.pkl` 保存已训练的 Model C，可直接用于预测及重新导出离散表；`udf/generated/model_c_table.h` 可直接供 Fluent UDF 编译。只加载可信来源的 pickle 文件。

## 从实验重新训练

需要在本地补齐 `data/raw/rheology/Q1.csv` 至 `Q7.csv`。这些路径被 `.gitignore` 排除。输入格式以 `src/data/rheology.py` 为准，单位和样品映射以 `configs/model_c.json` 为准。没有原始曲线时，不能独立重现训练和实验误差指标。

本仓库的图与指标是现有运行结果，不代表本次整理重新进行了训练或设备验证。模型配置含 PVDF 参数；用于 CMC 水体系前需核对适用性。
