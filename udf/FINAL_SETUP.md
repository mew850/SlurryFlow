# Model C + Mixture + CB-PBM 最终工程实现

本版实现固定分形维度的聚集、破碎、PSD结构反馈及Model C查表。
它是可标定的工程闭合模型，不代表核参数、分形维度或设备结果已实验验证。
旧文件保留；不要把旧黏度/Adjust/PBM核与本版混合挂接。

## 1. 文件和编译

Source Files只添加 `model_c_pbm_final.c`。
Header Files添加 `cb_final_config.h`、`cb_final_math.h`、`model_c_table.h`。
后三个文件与C文件放在同一个英文目录。项目中的表头原件在generated子目录。
不要把项目检查生成的ud_io1.h手工添加到Fluent Header Files；Fluent构建器会生成它。

使用双精度CPU Fluent；本机25.1真实头文件的3ddp_host和3ddp_node语法检查已通过。
未通过GUI对用户网格进行加载和瞬态求解验证；本包不提供跨版本预编译DLL。

## 2. 适用设置

- Mixture，恰好3相：胶液0、石墨1、CB2。
- 两个次相non-granular；无空气，无滑移。当前不覆盖自由表面、脱泡、逐次加料。
- PBM = homogeneous Discrete，选择CB，恰好15个Bin。
- CB真密度用于固体质量；PBM直径是固体体积等效直径dv。
- 不开启Growth/Nucleation/质量传递。
- 主相材料上的黏度是有效流变闭合；核对Fluent mixture laminar viscosity与UDM3，
  不能只因为UDF返回正确就认定动量方程使用的混合物黏度相同。
  若两者明显不同，应先审查相应力/混合物黏度处理，禁止直接用于定量力矩预测。
- 当前Model C表不以局部石墨浓度为独立输入。CMC水体系的参数适用性仍需核对。

## 3. 模型

dp为基本不可破碎单元直径，df固定，分形前因子取1。

da = dp * (dv/dp)^(3/df)

S = sum_i f_i * (dv_i/dp)^(9/df-3)

R = S/Sref

phi_input = alpha_CB * R^lambda

eta = ModelC(phi_input, gamma)

其中fi为CB内部固体体积分额；C_PB_DISCI原始值先检查再归一化。
代码接受原始和近似1或alpha_CB两种尺度，记录相对误差；偏差超过1%返回状态-4，
该单元临时回退R=1，必须排错，不能把错误状态算例用于结果分析。
PSD为NaN、明显负数、空分布或存储缺失时同样回退并标记。
UDM10保留原始和；UDM22记录归一化前误差。

聚集：Es*gamma*(da1+da2)^3/6 + EB*2*kB*T/(3*muB)*(da1+da2)*(1/da1+1/da2)。
与旧版相比，碰撞直径改成外包络da，因此聚集速度可能显著变化，需要重新标定。
默认muB使用当前Model C有效黏度；配置可改用实测胶液黏度，这是不同阻力闭合假设。

破碎：Cb*gamma_ref*(gamma/gamma_ref)^p * max((dv/dref)^q-1,0)。
本版破碎阈值按dv定义。子颗粒采用内置等体积二元分裂。
最小颗粒以下的分裂被禁止；最大粒径碰撞仍可能越界，因此必须检查边界Bin堆积。

本模型没有时间演化的df，没有额外UDS；记忆由PBM输运PSD提供。
Property每次从当前PSD计算，Adjust和Execute-at-end刷新诊断；没有按迭代次数累积的结构更新。

## 4. 默认参数和参考结构

配置文件cb_final_config.h：

| 参数 | 默认值 | 说明 |
|---|---:|---|
| CB_FINAL_DF | 2.3 | 演示分形维度，不等同于Model C拟合参数的验证结果 |
| CB_FINAL_DP | 1e-7 m | 演示基本单元直径，不能大于最小Bin |
| CB_FINAL_REF_DV | 8.483428982e-7 m | 参考单分散质量等效直径 |
| CB_FINAL_SREF_OVERRIDE | 0 | >0时使用用户提供的参考PSD结构系数 |
| CB_FINAL_FEEDBACK | 1 | 默认开启反馈，0为对照 |
| CB_FINAL_LAMBDA | 0.3 | 调试反馈强度，后续校准 |
| CB_FINAL_PROBE_GAMMA | 10 s^-1 | 固定测试剪切率 |
| CB_FINAL_ES / EB | 0.5 / 1 | 剪切/Brownian碰撞效率 |
| CB_FINAL_BREAK_C | 0.01 | 破碎强度 |
| CB_FINAL_BREAK_DREF | 1e-6 m | 质量等效破碎阈值 |
| CB_FINAL_BREAK_GREF | 100 s^-1 | 参考剪切率 |

默认Sref=7.04410931，对应参考外包络直径da=1.62622309e-6 m。
参考PSD不是每个单元随时间重置的量。本版使用固定参考，重启保持同一配置即可。
非均匀初始PSD下R不必等于1；不要按每个单元捕获不同参考而丢掉初始结构差异。

已有算例若dv初始分布不同，可以继续求解，但必须明确参考结构，不能把原Case当前PSD
自动当作流变标定样品PSD。使用实测PSD时先把测得外包络尺度转换为dv，并明确分布权重。

## 5. Bin设置和编号验证（必须）

输入原configs/fluent/pbm_cb_bins_15.txt，15组0.1--39.86微米，含义改为dv。
默认 `CB_FINAL_BIN0_LARGEST=1`，即Bin0最大、Bin14最小。
如果实际编号相反，改成0再编译。
执行On Demand的cb_final_check会打印本UDF假设的全部Bin编号、dv、da。
与Fluent的实际定义逐项核对。File模式若没有Print Bins按钮，使用PBM数密度输出
或单Bin初始化试验核对，不要只根据文件行序推断编号。

默认参考dv是升序第6个值；在默认降序映射中为Bin9。
单Bin9初始化后应得到UDM16=8.483428982e-7、UDM17=1.62622309e-6、UDM13=1。
务必同时用Fluent自带PSD/矩计算核对，单独用UDF自身输出不能验证其编号假设。
若更改Bin数或粒径，必须同步配置数组和Fluent；仅15组存储存在不代表粒径值匹配。

## 6. 迁移及挂接

1. 保存旧Case/Data。新建副本，移除旧Adjust、旧Execute-at-end及旧属性/核挂接。
2. 在Mixture层级分配24个UDM，建议全域分配，不使用zone-specific UDM allocation。
3. 以新库名libudf_cb_final编译本包唯一C源文件及三个自定义头文件。
4. 主相Materials/Viscosity挂cb_final_viscosity；不要挂在CB或石墨相。
5. Function Hooks/Adjust挂cb_final_adjust。
6. Function Hooks/Execute At End挂cb_final_end。
7. PBM Aggregation挂cb_final_aggregation，Factor=1。
8. PBM Breakage Frequency挂cb_final_breakage，Factor=1。
9. Breakage PDF=generalized，One Term，Daughters=2，qi0=1e20，Validate/Apply。
10. cb_final_check为按需诊断，不是初始化函数，不写体积分数或PSD。

所有DEFINE声明均为单行。旧版若同库编译虽函数名不同，也会增加误挂接风险。

## 7. 初始化与开算

推荐新验证从均匀CB=0.0152和实际石墨体积分数开始，主相为余量。
保持原已验证的立方体边界；可用0.01 m边长、顶壁1 m/s的剪切驱动示例。
PBM参考单分散：核对编号后把参考dv对应Bin设1，其余0，内部Bin和=1。
初始化后运行少量迭代/时间步使Data_Valid_P和梯度可用，再执行cb_final_check。
检查UDM4=1、UDM13=1（参考初始化），UDM9约1，UDM22很小。
UDM4=0是尚未有有效解，不能作为运行成功的标志。

从dt=1e-4--1e-3 s试算。分形碰撞尺寸比旧版大，不能照搬旧时间步。
根据源项刚性、每步Bin变化、CFL和内迭代收敛调整；dt减半检查是必须的。
先FEEDBACK=0跑对照，再FEEDBACK=1；两次必须从相同初始PSD/流场和相同物理时间起点开始。
调试lambda可从0.1增加到0.3，不要把递增lambda产生的瞬变解释为真实物理过程。

## 8. UDM和状态码

| UDM | 内容 |
|---:|---|
| 0/1 | CB/石墨体积分数 |
| 2 | 未截断等效查表浓度 |
| 3 | 当前局部剪切率下Model C黏度 |
| 4 | 状态码 |
| 5 | Mixture剪切率 |
| 6/7 | 浓度/剪切率越界，-1低于范围、1高于范围 |
| 8/9 | 主相体积分数/三相和 |
| 10 | 原始PBM变量和 |
| 11/12/13 | S/Sref/R |
| 14/15 | 原始/实际查表浓度 |
| 16/17 | 质量等效D32 / 外包络D32，m |
| 18 | 最大粒径组内部体积分额 |
| 19 | 浓度限幅标记0或1 |
| 20 | 固定10 s^-1测试剪切率下黏度，Pa.s |
| 21 | 根据CB质量和归一化PSD重构的颗粒数密度，1/m^3 |
| 22 | 原始Bin和的相对误差 |
| 23 | 最小粒径组内部体积分额 |

状态1=正常，2=几乎无CB，0=初始化/解未就绪；负值均须排错：
-1相数量/索引，-2 PBM存储缺失，-3 PSD无效或明显负值，-4原始Bin和不一致，-5相体积分数无效。
错误状态下临时回退原Model C，不表示计算可以作为有效结果使用。
UDM21是重构诊断，不是独立守恒证明；另需Fluent原生矩和CB质量积分。

## 9. 必做验证

- 单Bin参考恢复：R=1、eta恢复旧表；同时核对Fluent原生粒径。
- FEEDBACK=0/1同起点对比。FEEDBACK=0只关闭PSD到黏度，碰撞核仍采用本版分形映射。
- 关闭所有聚集/破碎：均匀PSD应保持，结构不应漂移。
- 仅聚集：全域总数下降；仅破碎：上阈值初始颗粒总数上升。
- 两种情况CB总质量、Fluent原生三阶矩基本守恒；最大/最小Bin无明显截断堆积。
- 监控UDM6、19：大量浓度越界说明表范围/反馈假设不足，不是结构饱和。
- dt减半、增加Bin数（同步配置）、网格加密，比较同物理时刻的PSD和力矩/壁面力。
- 低剪切恢复实验比较UDM20，避免把降低剪切率导致的瞬时增黏当作结构回弹。
- 最终与流变/PSD时间序列标定，不承诺给定参数一定出现回弹或动态平衡。

## 10. 检查和来源

本地check_final_compile.ps1使用安装的Fluent 2025 R1头文件和内置clang对host/node检查，
不是DLL链接或Fluent案例运行。tests/test_cb_final_math.c使用同一数学头文件测试
全部16281表节点、限幅、参考恢复、分形质量关系、等体积子颗粒和关闭反馈恢复。

官方接口/设置：https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/flu_ug_models_popbal.html
分形映射、参考归一化反馈及本版破碎频率属于项目自定义闭合，不是官方PBM默认模型。
