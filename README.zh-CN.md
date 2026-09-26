# GitHub上传说明

该仓库包含分析代码、参数字典、标准化分析表、输出结果、复现记录，以及面向Phytoparasitica稿件重新设计的主图和补充图。

## 推荐操作

1. 打开PowerShell并进入本文件夹。
2. 运行 `./CHECK_BEFORE_UPLOAD.ps1`。
3. 如果已经安装并登录GitHub CLI，运行：

```powershell
./UPLOAD_TO_GITHUB.ps1 -Visibility public
```

脚本会初始化Git仓库、创建首次提交、在GitHub建立仓库并推送。默认仓库名是`cimmyt-wheat-blast-cross-year-selection`。

如果没有安装GitHub CLI，也可以在GitHub网页新建空仓库，然后上传本文件夹中的全部内容。

## 发布前需要替换

- GitHub发布`v1.0.0`后，建议连接Zenodo生成永久DOI。
- 获得Zenodo DOI后，更新`CITATION.cff`和论文Data Availability Statement。

## 数据说明

仓库中的CIMMYT源文件及其衍生数据不受本仓库MIT许可证重新授权，其使用和再分发仍受原数据集条款约束。原数据链接列于`DATA_LICENSE.md`。
