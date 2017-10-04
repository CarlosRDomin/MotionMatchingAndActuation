function saveFigToFile(outputFileName, figHandle)
	if nargin<2 || isempty(figHandle), figHandle = gcf; end
	figFolderName = [fileparts(mfilename('fullpath')) '/../../figures/'];
	
	savefig(figHandle, [figFolderName outputFileName '.fig']);
	saveas(figHandle, [figFolderName outputFileName '.eps'], 'epsc');
end
