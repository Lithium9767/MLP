#!/usr/bin/env python3
"""Export a source-backed comparison of frozen discovery/validation profiles."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--clusters', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with args.clusters.open(encoding='utf-8-sig', newline='') as handle:
        rows = sorted(csv.DictReader(handle), key=lambda r: int(r['cluster']))
    if [int(r['cluster']) for r in rows] != list(range(12)):
        raise ValueError('Expected all fixed cluster IDs 0-11')
    profiles = [np.array([json.loads(r[f'{split}_reference_profile_1_66_and_unmapped']) for r in rows])
                for split in ('discovery', 'validation')]
    if any(x.shape != (12, 67) or not np.isfinite(x).all() or not np.allclose(x.sum(1), 1) for x in profiles):
        raise ValueError('Profiles must be 12 normalized 67-position distributions')
    tv = abs(profiles[0]-profiles[1]).sum(1)/2
    if not np.allclose(tv, [float(r['reference_profile_TV']) for r in rows]):
        raise ValueError('Displayed profiles and TV column disagree')
    plt.rcParams.update({'font.family': 'Microsoft YaHei', 'font.size': 11,
                         'axes.unicode_minus': False, 'figure.facecolor': 'white'})
    fig, axes = plt.subplots(1, 3, figsize=(15, 6.6), sharey=True,
                              gridspec_kw={'width_ratios': [4.8, 4.8, 1.7]})
    vmax = max(float(x.max()) for x in profiles)
    vmax = max(.1, np.ceil(vmax/.05)*.05)
    titles = ['发现集：1,202 条 / 12,175 窗口', '冻结迁移验证：519 条 / 5,276 窗口']
    for axis, values, title in zip(axes[:2], profiles, titles):
        drawing = axis.imshow(values, aspect='auto', cmap='Blues', vmin=0, vmax=vmax,
                              extent=(.5, 67.5, 11.5, -.5), interpolation='nearest')
        axis.set_title(title, pad=36, fontsize=12)
        axis.set_xticks([1, 11, 24, 35, 50, 62, 67], ['1', '11', '24', '35', '50', '62', '未映射'])
        axis.set_xlabel('7R1C 参考残基位置')
        for position in (10.5, 21.5, 23.5, 34.5, 36.5, 48.5, 49.5, 62.5, 66.5):
            axis.axvline(position, color='#777777', linewidth=.5, alpha=.5)
        for center, text in ((16, 'α1'), (29, 'β1'), (35.5, '转角'), (42, 'β2'), (56, 'α2')):
            axis.text(center, 1.025, text, transform=axis.get_xaxis_transform(), ha='center', fontsize=10)
    axes[0].set_yticks(range(12), [f'C{i}' for i in range(12)])
    axes[0].set_ylabel('冻结的簇编号')
    axes[2].barh(range(12), tv, height=.64, color='#CB8D3C')
    axes[2].set_xlim(0, .215)
    axes[2].set_title('位置分布差异 TV', pad=36, fontsize=12)
    axes[2].set_xlabel('0＝相同，1＝无重叠')
    axes[2].set_xticks([0, .1, .2])
    for i, value in enumerate(tv):
        label = '<0.001' if 0 < value < .001 else f'{value:.3f}'
        axes[2].text(value+.004, i, label, va='center', fontsize=10)
    for axis in axes:
        axis.spines[['top', 'right']].set_visible(False)
    fig.suptitle('K＝12：保留原模型后，区域对应关系能否在验证集复现',
                 x=.04, y=.98, ha='left', fontsize=16, fontweight='bold')
    fig.subplots_adjust(left=.06, right=.96, top=.81, bottom=.24, wspace=.23)
    color_axis = fig.add_axes([.07, .132, .46, .022])
    bar = fig.colorbar(drawing, cax=color_axis, orientation='horizontal', ticks=np.linspace(0, vmax, 4))
    bar.ax.xaxis.set_major_formatter(PercentFormatter(1))
    bar.set_label('每簇窗口在该参考位置的残基比例（全部 10 个残基为分母）', fontsize=10)
    fig.text(.62, .09, 'C10 的混合区域变化最大\n固定簇号比较；未重新聚类或匹配标签', fontsize=11)
    fig.text(.06, .035, '簇内按同源组 → 序列 → 窗口等权。此图检验位置分布复现，不能解释为功能预测准确率。', fontsize=10)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180, facecolor='white')
    plt.close(fig)
    record = {'source_csv_sha256': hashlib.sha256(args.clusters.read_bytes()).hexdigest(),
              'figure_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'matplotlib': matplotlib.__version__, 'shared_scale': [0, vmax],
              'figure_rows_match_source_TV': True, 'function_accuracy_displayed': False}
    args.output.with_suffix('.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
