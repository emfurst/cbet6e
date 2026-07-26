# Exporting notebooks to PDF

The file `vdw_report.tex.j2` formats the LaTeX PDF output to use a specific filename. It makes additional changes to the format of the code blocks.

To create the output, use the `nbeconvert` utility in `jupyter`:

```
% jupyter nbconvert van_der_waals_EOS.ipynb  --to pdf  --TemplateExporter.extra_template_paths=.  --template-file vdw_report.tex.j2 --no-prompt
```

