%% 设置导入选项并导入数据
opts = spreadsheetImportOptions("NumVariables", 2);

% 指定工作表和范围
opts.Sheet = "Datas for multi-period OPF";
opts.DataRange = "C2:D25";

% 指定列名称和类型
opts.VariableNames = ["Load", "wind"];
opts.VariableTypes = ["double", "double"];

% 导入数据
Wind = readtable("DATA.xlsx", opts, "UseExcel", false);
Wind=table2array(Wind);

%% 清除临时变量
clear opts