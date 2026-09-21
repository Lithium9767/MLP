# 蛋白语言模型表示

主责：D；复核：E。

主线采用ESM-2残基表示：

- 当前参考模型：`facebook/esm2_t12_35M_UR50D`；
- 残基向量：480维；
- 全长表示：去除特殊token后mean pooling；
- 局部表示：30 aa窗口、步长5；
- 缓存以sequence hash和模型版本命名，不提交Git。

模型或层改变时必须同步记录真实维度、权重版本、池化、设备和输出哈希。全长表示主要用于检查长度/物种等混杂；局部表示用于区域分析。

迁移后的入口为 `features/esm2_embed.py`。它默认只读取discovery primary记录，并拒绝在没有显式 `--allow-validation` 时打开validation：

```powershell
python features/esm2_embed.py `
  --metadata data/processed/gvpa_v1/metadata.csv `
  --split-manifest data/processed/gvpa_v1/split/split_manifest.csv `
  --output-dir embeddings/m3_discovery
```

可选运行依赖见 `requirements-m3.txt`。embedding和逐条manifest不提交Git。
