%[text] 内置anovan的升级版，将变量表作为分组输入，额外支持多重比较
%[text] ## 语法
%[text] ```matlabCodeExample
%[text] varargout=StatisticsAndMachineLearning.TabularAnovaN(Y,GroupTable,Name=Value);
%[text] %非平衡采样，指定一列因变量和自变量表
%[text] 
%[text] varargout=StatisticsAndMachineLearning.TabularAnovaN(YColumn,GroupTable,Name=Value);
%[text] %非平衡采样，指定因变量名和所有变量表
%[text] 
%[text] varargout=StatisticsAndMachineLearning.TabularAnovaN(OneWayTable,Name=Value);
%[text] %平衡采样，执行单因素方差分析
%[text] 
%[text] varargout=StatisticsAndMachineLearning.TabularAnovaN(NWayTensor,SampleDimension,Name=Value);
%[text] %平衡采样，执行多因素方差分析
%[text] ```
%[text] ## 输入参数
%[text] #### Y
%[text] 同内置anovan的第一个位置参数
%[text] #### GroupTable table
%[text] 分组变量，每列一个变量，每行对应一个Y
%[text] #### YColumn
%[text] GroupTable中要作为Y值的列名或列序号
%[text] #### OneWayTable table
%[text] 平衡采样的单因素数据表，第1维采样，第2维分组名称，允许组间单因素比较
%[text] #### NWayTensor
%[text] 平衡采样的多因素张量
%[text] #### SampleDimension(1,1)
%[text] 指定NWayTensor的哪个维度视为采样维度，其它维度视为因素
%[text] ### 名称值参数
%[text] #### Alpha(1,1)
%[text] 置信界限的显著性水平，指定为范围为0到1的标量值。
%[text] #### Continuous(1,:)
%[text] 连续预测因子的指示器，表示哪些分组变量应被视为连续预测因子，而不是分类预测因子，指定为变量名字符串向量。输入OneWayTable时，忽略此参数。
%[text] #### Display(1,1)logical=true
%[text] 显示方差分析表的指示器，当设为false时，只返回输出参数，而不将标准ANOVA表显示为图形。
%[text] #### Model='linear'
%[text] 模型的类型，指定为以下选项之一：
%[text] - 'linear'，默认的线性模型只计算N个主要效应的零假设的p值。
%[text] - 'interaction'，交互模型计算N个主效应和$C\_N^2${"editStyle":"visual"}双因素相互作用。
%[text] - 'full'，完整模型计算所有水平上N个主要效应和相互作用的零假设的p值。
%[text] - (1,1)uint8，对于k的整数值，(k≤N)为模型类型，计算到第k层的所有交互级别。例如，值3表示主效应加上两因素和三因素相互作用。k=1和k=2分别相当于线性和交互模型。值k=N相当于完整模型。
%[text] - tabular，每一列都是(:,1)logical的表，列名是变量名。为了更精确地控制ANOVA计算的主项和交互项，您可以指定一个表格，其中包含要包含在ANOVA模型中的每个主项或交互项的一行。每一行用N个true和false的定义一项，指示每个变量是否包含在项中。未列出的变量一律视为false。 \
%[text] 输入OneWayTable时，忽略此参数。
%[text] #### Nested tabular
%[text] 行名和列名都是变量名，但不必包含所有变量。每一列都是(:,1)logical，指示组变量之间的嵌套关系。例如，如果变量i嵌套在变量j中，则Nested{i,j}=true。不能在连续变量中指定嵌套。平衡采样不支持此参数。
%[text] #### Random(1,:)string
%[text] 随机变量的指示器，表示哪些分组变量是随机的，指定为变量名字符串向量。默认情况下，将所有分组变量视为固定的。如果交互项中的任何变量是随机的，则将交互项视为随机的。
%[text] #### SSType(1,1)
%[text] 平方和的类型，指定为以下选项之一：
%[text] - 1，类型Ⅰ平方和。通过将该项添加到已经包含在它之前列出的项的拟合中而得到的残差平方和的减少。
%[text] - 2，类型Ⅱ平方和。通过将该项添加到由不包含该项的所有其他项组成的模型中而得到的残差平方和的减少。
%[text] - 3，类型Ⅲ平方和。通过将该项添加到包含所有其他项的模型中而获得的残差平方和的减少，但其效果受制于通常的“sigma限制”，使模型可估计。
%[text] - 'h'，层次模型。类似于类型Ⅱ，但使用连续和分类因素来确定项的层次结构。 \
%[text] 任何项的平方和都是通过比较两个模型来确定的。对于包含主效应但不包含相互作用的模型，sstype的值仅影响非平衡数据的计算。
%[text] #### Comparison
%[text] 要额外计算多重比较P值的对组。
%[text] 如果输入是table，则Comparison也必须是table，每行一个比较对组，每列是一个变量，列名就是变量名（不允许使用PValue作为列名），列值是(:,2)，在第2维上排列对组中该变量的两个值。
%[text] 如果输入的是NWayTensor，则Comparison必须是(2,ndims(NWayTensor)-1,:)张量。第1维是要比较的两组，第2维按顺序排列各因素维度，第3维是不同的比较组。采样维度应当跳过，因此Comparison第2维长度应比NWayTensor维度数少一个。
%[text] ## 返回值
%[text] 如果指定了Comparison参数：
%[text] - 如果输入是table，则第一个返回值也是table，是在输入的Comparison基础上加一列PValue，表示每个比较对组的多重比较P值。后续返回值同内置anovan。
%[text] - 如果输入是张量，则第一个返回值是(:,1)，与输入Comparison的第3维一一对应，表示每个比较对组的P值。后续返回值同内置anovan。 \
%[text] 如果未指定Comparison参数，则返回值同内置anovan。
%[text] **See also** [anovan](<matlab:doc anovan>) [StatisticsAndMachineLearning.TabularAnova1](<matlab:doc StatisticsAndMachineLearning.TabularAnova1>)
function varargout = TabularAnovaN(varargin)

[mode,positionalCount] = iDetectSyntax(varargin);

nv = varargin(positionalCount+1:end);
[AnovanOptions,MultcompareOptions,provided] = iParseNameValue(nv);

switch mode
	case "GroupTable"
		Y = varargin{1};
		GroupTable = varargin{2};

		% 允许用列序号或列名指定因变量
		if isnumeric(Y) && isreal(Y) && isscalar(Y)
			try
				YColumn = GroupTable{:,Y};
				GroupTable(:,Y) = [];
				Y = YColumn;
			catch
			end
		elseif isstring(Y) || ischar(Y)
			YColumn = GroupTable.(Y);
			GroupTable.(Y) = [];
			Y = YColumn;
		end

		VarNames = GroupTable.Properties.VariableNames;
		NumVariables = width(GroupTable);

		AnovanOptions = iNormalizeAnovanOptions(AnovanOptions,VarNames,NumVariables);
		AnovanOptionsCell = [lower(fieldnames(AnovanOptions)),struct2cell(AnovanOptions)]';

		Groups = cell(1,NumVariables);
		for V = 1:NumVariables
			Groups{V} = GroupTable{:,V};
			if isinteger(Groups{V})
				% anovan 不支持整数类型
				Groups{V} = categorical(Groups{V});
			end
		end
		if isduration(Y)
			Y = seconds(Y);
		end

		[varargout{1:clip(nargout,3,4)}] = anovan(Y,Groups,AnovanOptionsCell{:},varnames=VarNames);
		if provided.Comparison
			Comparison = MultcompareOptions.Comparison;
			[Comparison,ComparisonPValue] = iComparisonFromAnovanStats(varargout{3},Comparison,AnovanOptions,provided);
			% table 输入：返回加 PValue 的 table
			Comparison.PValue = ComparisonPValue;
			varargout = [{Comparison},varargout];
		end

	case "OneWayTable"
		OneWayTable = varargin{1};
		X = OneWayTable{:,:};
		if isduration(X)
			X = seconds(X);
		end

		% OneWayTable 语法：使用 anova1（连续变量等 anovan 选项不适用）
		alphaProvided = provided.Alpha;
		displayProvided = provided.Display;
		alphaVal = 0.05;
		if alphaProvided && ~isempty(AnovanOptions.Alpha)
			alphaVal = AnovanOptions.Alpha;
		end

		if displayProvided
			displayStr = 'on';
			if ~AnovanOptions.Display
				displayStr = 'off';
			end
			[p,tbl,stats] = anova1(X,OneWayTable.Properties.VariableNames,alphaVal,displayStr);
		else
			if alphaProvided
				[p,tbl,stats] = anova1(X,OneWayTable.Properties.VariableNames,alphaVal);
			else
				[p,tbl,stats] = anova1(X,OneWayTable.Properties.VariableNames);
			end
		end

		outCount = clip(nargout,3,4);
		varargout(1:outCount) = {p,tbl,stats,[]};

		if provided.Comparison
			Comparison = MultcompareOptions.Comparison;
			[Comparison,ComparisonPValue] = iComparisonFromAnova1Stats(stats,Comparison,AnovanOptions,provided,alphaVal);
			Comparison.PValue = ComparisonPValue;
			varargout = [{Comparison},varargout];
		end

	case "NWayTensor"
		NWayTensor = varargin{1};
		SampleDimension = varargin{2};

		if ~isnumeric(SampleDimension) || ~isscalar(SampleDimension) || SampleDimension < 1 || SampleDimension ~= fix(SampleDimension)
			StatisticsAndMachineLearning.Exception.SampleDimension_invalid.Throw;
		end

		nd = ndims(NWayTensor);
		if SampleDimension > nd
			StatisticsAndMachineLearning.Exception.SampleDimension_exceeds_ndims.Throw;
		end

		Y = NWayTensor;
		if isduration(Y)
			Y = seconds(Y);
		end
		Y = Y(:);

		sz = size(NWayTensor);
		idx = (1:numel(NWayTensor))';
		subs = cell(1,nd);
		[subs{:}] = ind2sub(sz,idx);

		factorDims = setdiff(1:nd,SampleDimension,'stable');
		numFactors = numel(factorDims);
		Groups = cell(1,numFactors);
		for k = 1:numFactors
			Groups{k} = subs{factorDims(k)};
		end
		VarNames = cellstr("Dim" + string(factorDims));

		if provided.Nested
			StatisticsAndMachineLearning.Exception.Nested_unsupported_balanced.Throw;
		end

		AnovanOptions = iNormalizeAnovanOptions(AnovanOptions,VarNames,numFactors);
		AnovanOptionsCell = [lower(fieldnames(AnovanOptions)),struct2cell(AnovanOptions)]';

		[varargout{1:clip(nargout,3,4)}] = anovan(Y,Groups,AnovanOptionsCell{:},varnames=VarNames);
		if provided.Comparison
			ComparisonTensor = MultcompareOptions.Comparison;
			ComparisonTable = iTensorComparisonToTable(ComparisonTensor,VarNames);
			[~,PValue] = iComparisonFromAnovanStats(varargout{3},ComparisonTable,AnovanOptions,provided);
			% 张量输入：第一个返回值为 PValue 向量
			varargout = [{PValue},varargout];
		end
end

end

function [mode,positionalCount] = iDetectSyntax(args)
if isempty(args)
	StatisticsAndMachineLearning.Exception.Not_enough_input_arguments.Throw;
end

if istable(args{1}) && (isscalar(args) || ~istable(args{2}))
	mode = "OneWayTable";
	positionalCount = 1;
	return
end

if numel(args) >= 2 && istable(args{2})
	mode = "GroupTable";
	positionalCount = 2;
	return
end

if numel(args) >= 2 && isnumeric(args{2}) && isscalar(args{2})
	mode = "NWayTensor";
	positionalCount = 2;
	return
end

	StatisticsAndMachineLearning.Exception.Input_syntax_unsupported.Throw;
end

function [AnovanOptions,MultcompareOptions,provided] = iParseNameValue(nv)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = false;

addParameter(ip,'Alpha',[]);
addParameter(ip,'Continuous',[]);
addParameter(ip,'Display',[]);
addParameter(ip,'Model',[]);
addParameter(ip,'Nested',[]);
addParameter(ip,'Random',[]);
addParameter(ip,'SSType',[]);
addParameter(ip,'Comparison',[]);

parse(ip,nv{:});
r = ip.Results;

provided = struct;
provided.Alpha = ~isempty(r.Alpha);
provided.Continuous = ~isempty(r.Continuous);
provided.Display = ~isempty(r.Display);
provided.Model = ~isempty(r.Model);
provided.Nested = ~isempty(r.Nested);
provided.Random = ~isempty(r.Random);
provided.SSType = ~isempty(r.SSType);
provided.Comparison = ~isempty(r.Comparison);

AnovanOptions = struct;
if provided.Alpha, AnovanOptions.Alpha = r.Alpha; end
if provided.Continuous, AnovanOptions.Continuous = r.Continuous; end
if provided.Display, AnovanOptions.Display = r.Display; end
if provided.Model, AnovanOptions.Model = r.Model; end
if provided.Nested, AnovanOptions.Nested = r.Nested; end
if provided.Random, AnovanOptions.Random = r.Random; end
if provided.SSType, AnovanOptions.SSType = r.SSType; end

MultcompareOptions = struct;
if provided.Comparison, MultcompareOptions.Comparison = r.Comparison; end
end

function AnovanOptions = iNormalizeAnovanOptions(AnovanOptions,VarNames,NumVariables)
% Continuous
if isfield(AnovanOptions,'Continuous')
	if isstring(AnovanOptions.Continuous) || ischar(AnovanOptions.Continuous)
		[~,AnovanOptions.Continuous] = ismember(AnovanOptions.Continuous,VarNames);
	end
end

% Display
if isfield(AnovanOptions,'Display')
	if AnovanOptions.Display
		AnovanOptions.display = 'on';
	else
		AnovanOptions.display = 'off';
	end
end

% Model (tabular)
if isfield(AnovanOptions,'Model') && istabular(AnovanOptions.Model)
	modelTable = AnovanOptions.Model;
	modelTable{:,setdiff(VarNames,modelTable.Properties.VariableNames)} = false;
	AnovanOptions.Model = double(modelTable{:,VarNames});
	if any(~ismember(AnovanOptions.Model,[0,1]),'all')
		StatisticsAndMachineLearning.Exception.AnovaN_model_invalid.Throw;
	end
end

% Nested (tabular)
if isfield(AnovanOptions,'Nested') && istabular(AnovanOptions.Nested)
	Nested = AnovanOptions.Nested;
	AnovanOptions.Nested = false(NumVariables);
	[~,RowIndex] = ismember(Nested.Properties.RowNames,VarNames);
	[~,ColumnIndex] = ismember(Nested.Properties.VariableNames,VarNames);
	AnovanOptions.Nested(RowIndex,ColumnIndex) = Nested{:,:};
end

% Random
if isfield(AnovanOptions,'Random')
	if isstring(AnovanOptions.Random) || ischar(AnovanOptions.Random)
		[~,AnovanOptions.Random] = ismember(AnovanOptions.Random,VarNames);
	end
end
end

function [Comparison,ComparisonPValue] = iComparisonFromAnovanStats(stats,Comparison,AnovanOptions,provided)
if ~istable(Comparison)
	StatisticsAndMachineLearning.Exception.Comparison_must_be_table.Throw;
end

if any(strcmpi(Comparison.Properties.VariableNames,'PValue'))
	StatisticsAndMachineLearning.Exception.Comparison_contains_PValue.Throw;
end

AnovanOptionsCell = {};
if isfield(AnovanOptions,'display')
	AnovanOptionsCell = [AnovanOptionsCell,{'display',AnovanOptions.display}];
end
if provided.Alpha && isfield(AnovanOptions,'Alpha')
	AnovanOptionsCell = [AnovanOptionsCell,{'Alpha',AnovanOptions.Alpha}];
end

VariableNames = Comparison.Properties.VariableNames;
[~,Dimensions] = ismember(VariableNames,stats.varnames);
if any(Dimensions == 0)
	StatisticsAndMachineLearning.Exception.Comparison_varnames_mismatch.Throw;
end

[ComparisonMatrix,~,~,GNames] = multcompare(stats,AnovanOptionsCell{:},Dimension=Dimensions);

Groups = split(GNames,"="|",");
Groups = cell2table(Groups(:,2:2:end),VariableNames=Groups(1,1:2:end));

for C = string(VariableNames)
	Groups.(C) = iCoerceGroupColumnType(Groups.(C),Comparison.(C));
end

[ComparisonA,ComparisonB] = deal(table);
for C = 1:width(Comparison)
	ComparisonA.(VariableNames{C}) = Comparison{:,C}(:,1);
	ComparisonB.(VariableNames{C}) = Comparison{:,C}(:,2);
end

[Logical,ComparisonAIndex] = ismember(ComparisonA,Groups);
if ~all(Logical)
	StatisticsAndMachineLearning.Exception.Nonexistent_comparison_groups.Throw(ComparisonA(~Logical,:));
end
[Logical,ComparisonBIndex] = ismember(ComparisonB,Groups);
if ~all(Logical)
	StatisticsAndMachineLearning.Exception.Nonexistent_comparison_groups.Throw(ComparisonB(~Logical,:));
end

[~,rowIndex] = ismember(sort([ComparisonAIndex,ComparisonBIndex],2),ComparisonMatrix(:,1:2),'rows');
if any(rowIndex == 0)
	StatisticsAndMachineLearning.Exception.Nonexistent_comparison_groups.Throw(Comparison(rowIndex==0,:));
end
ComparisonPValue = ComparisonMatrix(rowIndex,6);
end

function [Comparison,ComparisonPValue] = iComparisonFromAnova1Stats(stats,Comparison,AnovanOptions,provided,alphaVal)
if ~istable(Comparison)
	StatisticsAndMachineLearning.Exception.OneWay_Comparison_must_be_table.Throw;
end

if any(strcmpi(Comparison.Properties.VariableNames,'PValue'))
	StatisticsAndMachineLearning.Exception.Comparison_contains_PValue.Throw;
end

if width(Comparison) ~= 1
	StatisticsAndMachineLearning.Exception.OneWay_Comparison_must_single_variable.Throw;
end

if isfield(AnovanOptions,'display')
	displayStr = AnovanOptions.display;
else
	displayStr = 'on';
end

mcArgs = {'Display',displayStr};
if provided.Alpha
	mcArgs = [mcArgs,{'Alpha',alphaVal}];
end

[ComparisonMatrix,~,~,GNames] = multcompare(stats,mcArgs{:});

groupVar = Comparison.Properties.VariableNames{1};
groupPairs = Comparison{:,1};
if size(groupPairs,2) ~= 2
	StatisticsAndMachineLearning.Exception.OneWay_Comparison_column_invalid.Throw;
end

% 允许用组名（字符串/字符/分类）或用数值索引指定对组
gn = string(GNames);
groupIndexA = iMapOneWayGroupToIndex(groupPairs(:,1),gn);
groupIndexB = iMapOneWayGroupToIndex(groupPairs(:,2),gn);

[~,rowIndex] = ismember(sort([groupIndexA,groupIndexB],2),ComparisonMatrix(:,1:2),'rows');
if any(rowIndex == 0)
	StatisticsAndMachineLearning.Exception.Nonexistent_comparison_groups.Throw(Comparison(rowIndex==0,:));
end

ComparisonPValue = ComparisonMatrix(rowIndex,6);

% 对 OneWayTable：把变量名统一命名为 Group，便于文档一致
if ~strcmp(groupVar,'Group')
	Comparison = renamevars(Comparison,groupVar,'Group');
end
end

function idx = iMapOneWayGroupToIndex(values,groupNames)
if isnumeric(values)
	idx = values;
	return
end

v = string(values);
[tf,idx] = ismember(v,groupNames);
if ~all(tf)
	StatisticsAndMachineLearning.Exception.Nonexistent_comparison_groups.Throw(table(v(~tf),'VariableNames',{'Group'}));
end
end

function out = iCoerceGroupColumnType(groupTokens,comparisonColumn)
% groupTokens: string array tokens from multcompare gnames
% comparisonColumn: (:,2) array in Comparison table

% comparisonColumn may be cell, categorical, logical, numeric, string
if islogical(comparisonColumn)
	out = (groupTokens=="1") | strcmpi(groupTokens,"true");
	return
end

if isnumeric(comparisonColumn)
	out = str2double(groupTokens);
	out = cast(out,class(comparisonColumn));
	return
end

if iscategorical(comparisonColumn)
	out = categorical(groupTokens);
	return
end

% default: string
out = string(groupTokens);
end

function ComparisonTable = iTensorComparisonToTable(ComparisonTensor,VarNames)
if ~(isnumeric(ComparisonTensor) || islogical(ComparisonTensor) || isstring(ComparisonTensor) || ischar(ComparisonTensor) || iscategorical(ComparisonTensor))
	StatisticsAndMachineLearning.Exception.Comparison_tensor_invalid_type.Throw;
end

if ndims(ComparisonTensor) ~= 3
	StatisticsAndMachineLearning.Exception.Comparison_tensor_invalid_ndims.Throw;
end

if size(ComparisonTensor,1) ~= 2
	StatisticsAndMachineLearning.Exception.Comparison_tensor_first_dim_invalid.Throw;
end

numFactors = numel(VarNames);
if size(ComparisonTensor,2) ~= numFactors
	StatisticsAndMachineLearning.Exception.Comparison_tensor_second_dim_mismatch.Throw;
end

K = size(ComparisonTensor,3);
ComparisonTable = table;
for k = 1:numFactors
	ComparisonTable.(VarNames{k}) = squeeze(permute(ComparisonTensor(:,k,:),[3 1 2]));
end

% 保证是 K 行
if height(ComparisonTable) ~= K
	ComparisonTable = ComparisonTable(1:K,:);
end
end

%[appendix]{"version":"1.0"}
%---
