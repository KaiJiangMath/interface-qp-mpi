%% plot morphologies in real space;
function plotPhase(model_type, bulk_phase1, bulk_phase2, PARA_flag, plot_flag, sflag)

	%clear;
	clc; close all;
	set(0, 'DefaultFigureVisible', 'off');
	if ( nargin == 5 )
		sflag = 'any';
	end
	% example: plotPhase('LB', 'BCC', 'BCC', '1', 'sys');
	%model_type = 'LB';
	%bulk_phase1 = 'DG';
	%bulk_phase2 = 'DG';
	%PARA_flag = '0';
	%plot_flag = 'sys';
	if ~( strcmp(plot_flag, 'bulk') || strcmp(plot_flag, 'sys') ||...
		strcmp(plot_flag, 'real') || strcmp(plot_flag, 'spec') ||...
	   	strcmp(plot_flag, 'all') )
		fprintf("\t       plot_flag only 'bulk', 'sys', 'real', 'spec', or 'all'\n");
		return;
	end

	fprintf('\t ***** model type: %s ***** \n', model_type);
	fprintf('\t ***** left bulk value: %s ***** \n', bulk_phase1);
	fprintf('\t ***** right bulk value: %s ***** \n', bulk_phase2);
	fprintf('\t ***** parameter case: %s ***** \n', PARA_flag);
	fprintf('\n\n');

	rsltDir = sprintf('%s_%s_%s_%s', model_type, bulk_phase1, bulk_phase2, PARA_flag);
	if exist(['./fig/', rsltDir]) == 0
		mkdir(['./fig/', rsltDir]);
	end
	if exist(['./vtk/', rsltDir]) == 0
		mkdir(['./vtk/', rsltDir]);
	end

	%% plot bulk phases;
	if ( strcmp(plot_flag, 'bulk') || strcmp(plot_flag, 'all') )
		optIter	  = -1;
		bulkStart = -1;
		bulkEnd	  = -2;
		sbulkparam1 = fn_bulk_main ( rsltDir, optIter, bulkStart, bulkEnd, 1 );
		sbulkparam2 = fn_bulk_main ( rsltDir, optIter, bulkStart, bulkEnd, 2 );
	end

	%% connect the initial interface structure;
	initDist = 0.0;
%	if ( ~isempty(sbulkparam1.comDensity) && ~isempty(sbulkparam2.comDensity) )
%		fn_disp_bulk_connect_density ( sbulkparam1, sbulkparam2, initDist );
%	end

	%% plot the interface structure;
	if ( strcmp(plot_flag, 'real') || strcmp(plot_flag, 'sys') || strcmp(plot_flag, 'all') )
		fprintf('\t ***** interface system ***** \n');
		ssysparam = struct();
		ssysparam.rsltDir = rsltDir;
		ssysparam = fn_obt_system_param ( ssysparam );

		if strcmp(sflag, 'recover')
			if exist(['./fig/', rsltDir, '/recover']) == 0
				mkdir(['./fig/', rsltDir, '/recover']);
			end
			step = -2;
			ssysparam = fn_disp_system_density ( ssysparam, step, -1e4, 'read' );
			fn_vtk_system_density ( ssysparam, step );
		else
			step = 0;
			ssysparam = fn_disp_system_density ( ssysparam, step, -1e4, 'read' );
			fn_vtk_system_density ( ssysparam, step );
			ssysparam0 = ssysparam;

			step = -1;
			ssysparam = fn_disp_system_density ( ssysparam, step, -1e4, 'read' );
			fn_vtk_system_density ( ssysparam, step );
			ssysparam1 = ssysparam;

			%% plot after screening;
%			ssysparam = fn_disp_system_density ( ssysparam, step, 0.6, 'read' );

			%% save screened data;
%			ssysparam = fn_save_system_density ( ssysparam, step, 0.6 );
%			ssysparam = fn_save_system_density ( ssysparam, step, 0.4 );
%			ssysparam = fn_save_system_density ( ssysparam, step, 0.2 );

			%% compute the width of grain boundary;
			step = -10;
			ssysparam.density.val = ssysparam1.density.val - ssysparam0.density.val;
			ssysparam = fn_disp_system_density ( ssysparam, step, -1e4, 'noread' );
			fn_vtk_system_density ( ssysparam, step );

			%% the width of grain boundary;
			tol = 1e-1;
			len = size(ssysparam.density.val, 1);
			maxerr = zeros(1, len);
			for ej = 1:1:len
				maxerr(ej) = max(abs(ssysparam.density.val(ej,:)));
			end
			ind = find(maxerr > tol);
			gbreg = zeros(1, 2);
			if ~isempty(ind)
				gbreg(1) = ssysparam.density.xx(min(ind));
				gbreg(2) = ssysparam.density.xx(max(ind));
			end
			step = -11;
			ssysparam.density.val = ssysparam1.density.val;
			ssysparam = fn_disp_system_density ( ssysparam, step, -1e4, 'noread', gbreg );
			fn_vtk_system_density ( ssysparam, step );
		end
	end

	%% analyze interface spectra;
	if ( strcmp(plot_flag, 'spec') || strcmp(plot_flag, 'sys') || strcmp(plot_flag, 'all') )
		fprintf('\t ***** interface system ***** \n');
		ssysparam = struct();
		ssysparam.rsltDir = rsltDir;
		ssysparam = fn_obt_system_param ( ssysparam );

		%% obtain field information;
		TOL = 1.0e-2;
		step = -1;
		ssysparam = fn_obt_system_field ( ssysparam, step );

		if strcmp(sflag, 'recover')
			%% obtain 'rhoJCplx' after recovering;
			fn_obt_system_recover_density ( ssysparam, step, TOL );
		end

		%% sort spectra according to the intensity;
		if ( strcmp(sflag, 'sort') || strcmp(sflag, 'any') )
			[specIndComm, specIndAll] = fn_sort_system_field ( ssysparam, step, TOL );
		end

		%% analyze interface spectra on slice of 'xslice';
		if ( strcmp(sflag, 'xslice') || strcmp(sflag, 'any') )
			if exist(['./fig/', rsltDir, '/xslice']) == 0
				mkdir(['./fig/', rsltDir, '/xslice']);
			end
			xmod = 5; % please odd number;
			num = ssysparam.LGL_num + 1;
			fprintf("\t LGL_num = %d.\n", num);
			xslice_tot = round( linspace(1, num, xmod) );
%			xslice_tot = 10:10:250;
			if strcmp(PARA_flag, 'GBx50')
				xslice_tot = [xslice_tot, 10:10:250, 251:1:256];
			end
			for i = 1:1:length(xslice_tot)
				xslice = xslice_tot(i);
				fn_disp_system_field_xslice ( ssysparam, step, xslice, TOL )
			end
		end

		%% plot the curves of interface field against x;
		if ( strcmp(sflag, 'xcurve') || strcmp(sflag, 'any') )
			if exist(['./fig/', rsltDir, '/xcurve']) == 0
				mkdir(['./fig/', rsltDir, '/xcurve']);
			end
			[specIndComm, specIndAll] = fn_sort_system_field ( ssysparam, step, TOL );
			specIndTot = specIndAll;
%			specIndTot = [2, 33];
			for i = 1:1:length(specIndTot)
				specInd = specIndTot(i);
				fn_disp_system_field_xcurve ( ssysparam, step, specInd );
			end
		end
	end
end


%% main code of bulk phase;
function sbulkparam = fn_bulk_main ( rsltDir, optIter, bulkStart, bulkEnd, sflag )

	%% define structure body;
	switch sflag
	case 1
		fprintf('\t ***** parameters of left bulk phase ***** \n');
	case 2
		fprintf('\t ***** parameters of right bulk phase ***** \n');
	end
	sbulkparam = struct();
	sbulkparam.sflag   = sflag;
	sbulkparam.rsltDir = rsltDir;

	%% obtain parameters;
	sbulkparam = fn_obt_bulk_param( sbulkparam );
	if ( bulkStart < 0 && bulkEnd < 0 )
		save_step = -1;
	else
		save_step = sbulkparam.save_step;
	end

	%% start: save_step * bulkStart;
	%% end:   save_step * bulkEnd;
	for plotIter = bulkStart:save_step:bulkEnd
		%% display density;
		sbulkparam = fn_disp_bulk_density ( sbulkparam, optIter, plotIter );
%		return;
		%% display Fourier spectra;
		sbulkparam = fn_disp_bulk_field	  ( sbulkparam, optIter, plotIter );
	end

	%% plot bulk phases after projecting;
	sbulkparam = fn_disp_bulk_proj_density	 ( sbulkparam );

	%% plot bulk phases in the common space;
	sbulkparam = fn_disp_bulk_common_density ( sbulkparam );
	fprintf('\n');
end


%% obtain bulk parameters;
function sbulkparam = fn_obt_bulk_param ( sbulkparam )

	fname = sprintf('./result/%s/parameter_init.dat', sbulkparam.rsltDir);
%	fname = sprintf('./result/%s/parameter_opt.dat', sbulkparam.rsltDir);
	sflag = sbulkparam.sflag;

	data = textread(fname, '%s', 'delimiter', '\n');
	for i = 1:1:length(data)
		dataTmp = split(data{i});
		keyword = sprintf('bulk%d_dimPhy', sflag);
		if ( strcmp(dataTmp{1}, keyword) )
			dimPhy = str2num( dataTmp{3} );
		end
		keyword = sprintf('bulk%d_dimCpt', sflag);
		if ( strcmp(dataTmp{1}, keyword) )
			dimCpt = str2num( dataTmp{3} );
		end
		keyword = sprintf('bulk%d_rcpBox', sflag); % default reciprocal box;
		if ( strcmp(dataTmp{1}, keyword) )
			for j = 1:1:dimCpt
				buff = split( data{i+j} );
				for k = 1:1:length(buff)
					rcpBox(j,k) = str2num(buff{k});
				end
			end
		end
		keyword = sprintf('bulk%d_projMat', sflag);
		if ( strcmp(dataTmp{1}, keyword) )
			for j = 1:1:dimPhy
				buff = split( data{i+j} );
				for k = 1:1:length(buff)
					projMat(j,k) = str2num(buff{k});
				end
			end
		end
		keyword = sprintf('bulk%d_rotateMat', sflag);
		if ( strcmp(dataTmp{1}, keyword) )
			for j = 1:1:dimPhy
				buff = split( data{i+j} );
				for k = 1:1:length(buff)
					rotateMat(j,k) = str2num(buff{k});
				end
			end
		end
		keyword = sprintf('bulk%d_translVec', sflag);
		if ( strcmp(dataTmp{1}, keyword) )
			for j = 1:1:dimPhy
				buff = split( data{i+j} );
				for k = 1:1:length(buff)
					translVec(j,k) = str2num(buff{k});
				end
			end
		end
		keyword = sprintf('bulk%d_save_step', sflag);
		if ( strcmp(dataTmp{1}, keyword) )
			save_step = str2num( dataTmp{3} );
		end
		keyword = 'nprocs';
		if ( strcmp(dataTmp{1}, keyword) )
			nprocs = str2num( dataTmp{3} );
		end
	end

	%% search the optimized reciprocal box to replace 'bulk_rcpBox';
	for i = 1:1:length(data)
		dataTmp = split(data{i});
		keyword = sprintf('bulk%d_opt_rcpBox', sflag);
		if ( strcmp(dataTmp{1}, keyword) )
			for j = 1:1:dimCpt
				buff = split( data{i+j} );
				for k = 1:1:length(buff)
					rcpBox(j,k) = str2num(buff{k});
				end
			end
		end
	end
		
	%% save to the structure body 'sbulkparam';
	sbulkparam.dimPhy	 = dimPhy;
	sbulkparam.dimCpt	 = dimCpt;
	sbulkparam.rcpBox    = rcpBox;
	sbulkparam.projMat	 = projMat;
	sbulkparam.rotateMat = rotateMat;
	sbulkparam.translVec = translVec;
	sbulkparam.save_step = save_step;
	sbulkparam.nprocs	 = nprocs;
end


%% display density;
function sbulkparam = fn_disp_bulk_density ( sbulkparam, optIter, plotIter )

	%% read data;
	fname = sprintf('./result/%s/bulk%d_opt%d_density%d.dat',...
			sbulkparam.rsltDir, sbulkparam.sflag, optIter, plotIter);
	if ( exist(fname) == 0 )
		fprintf('\t       %s not exist.\n', fname);
		sbulkparam.density = [];
		return;
	end
	sbulkparam.density = fn_myfile_read ( fname, sbulkparam.dimPhy );

	%% dimension;
	dimPhy = sbulkparam.dimPhy;
	fprintf('\t ====> plot density: dimPhy = %d, plotIter = %d \n',...
		dimPhy, plotIter);

	set(0, 'DefaultFigureVisible', 'on');
	%% set configuration of the figure
	figTitle = sprintf('%s:density%d:%d', sbulkparam.rsltDir, sbulkparam.sflag, plotIter);
	Pos = [200, 250, 800, 600];
	figure('Name', figTitle, 'NumberTitle', 'off', 'Position', Pos);
	set(gca, 'LineWidth', 2, 'FontName', 'Times New Roman', 'FontSize', 20);
	hold on;

	%% different plotting cases;
	if ( dimPhy == 1 || dimPhy == 2 )
		%% spatial grid;
		x		= sbulkparam.density.xx;
		density = sbulkparam.density.val;
		ncpt	= size(density);
		switch dimPhy
		case 1
			y = x;
			density = repmat(density, ncpt(end:-1:1)); % expand along y-direction;
		case 2
			y = sbulkparam.density.yy;
		end
		[yy, xx] = meshgrid(y, x);

		%% plotting
		pcolor(xx, yy, density);
		shading interp;
		box on;
		axis tight;

		%% create colorbar
		colormap('jet');
%		colormap('gray');
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(y), max(y)]);
	elseif ( dimPhy == 3 )
		%% screen;
		x = sbulkparam.density.xx;
		y = sbulkparam.density.yy;
		z = sbulkparam.density.zz;
		density = sbulkparam.density.val;
		density = real(density);
		minu = min(density(:));
		maxu = max(density(:));
		isoA = minu + 0.7*(maxu-minu);
%		isoA = 1;

		%% spatial grid;
		[yy, xx, zz] = meshgrid(y, x, z);
		tol  = 1.0 * [-1, 1, -1, 1, -1, 1];
		axBC = [min(x), max(x), min(y), max(y), min(z), max(z)] + tol;
%		fn_disp_bulk_lattice ( sbulkparam, 2.0 * axBC, plotIter );

		%% plot and modify;
		fprintf("isoA is %.4e\n", isoA);
		alpA = 0.8;
		patch(isosurface(xx, yy, zz, density, isoA), 'facecolor', [30, 144, 255]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'none');
		patch(isocaps(xx, yy, zz, density, isoA, 'enclose'), 'facecolor', [0, 206, 209]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'interp'); % flat, interp
		daspect([1, 1, 1]);
		camup([1, 0, 0]);
		campos([25, -55, 5]);
		camlight;
		lighting phong;
		box on;
		axis on;
		axis(axBC);
%		view([112, 35]);
%		view([0,90]);
		view(3);

		%% create colorbar
		colormap('jet');
%		colormap('gray');
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		zlabel('$Z$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);

%		box off;
%		axis off;
%		colorbar off;
	end

	%% saving
	%set(gca, 'LooseInset', [0,0,0,0]);
	figName = sprintf('./fig/%s/bulk%d_density%d.png', sbulkparam.rsltDir,...
					sbulkparam.sflag, plotIter);
	saveas(gca, figName);
end


%% plot field of bulk phase (Fourier spectra);
function sbulkparam = fn_disp_bulk_field ( sbulkparam, optIter, plotIter )

	%% parameters;
	dimPhy	  = sbulkparam.dimPhy;
	dimCpt	  = sbulkparam.dimCpt;
	rotateMat = sbulkparam.rotateMat;
	projMat	  = sbulkparam.projMat;
	rcpBox	  = sbulkparam.rcpBox;
	nprocs	  = sbulkparam.nprocs;
	rotateProjBoxMat = rotateMat.' * projMat * rcpBox; % R'PB;

	%% generate spectra;
	spectra = [];
	FourCoeff = [];

	%% read data;
	for myrank = 0:1:nprocs-1
		fname = sprintf('./result/%s/rank%d/bulk%d_opt%d_field%d.dat',...
				sbulkparam.rsltDir, myrank, sbulkparam.sflag, optIter, plotIter);
		if ( exist(fname) == 0 )
			fprintf('\t       %s not exist.\n', fname);
		end
		field = load(fname);
		[nr, nc] = size(field);

		%% multiply R'PB;
		for i = 1:1:nr
			val = norm(field(i,dimCpt+1:1:end));
			if ( val > 1.0e-10 )
				coordinate = field(i,1:1:dimCpt) * rotateProjBoxMat.';
				spectra(end+1,:) = coordinate;
				FourCoeff(end+1) = val;
			end
		end
	end

	%% sort descendly and save into the structure body;
	[val, ind]  = sort(FourCoeff, 'descend');
	sbulkparam.field = struct('spectra', spectra(ind,:), 'val', val);


	dimPhy = sbulkparam.dimPhy;
	fprintf('\t ====> plot field:  dimPhy = %d, plotIter = %d \n', dimPhy, plotIter);

	%% set configuration of the figure
	Pos = [200, 250, 800, 600];
	figTitle = sprintf('%s:field%d', sbulkparam.rsltDir, sbulkparam.sflag);
	figure('Name', figTitle, 'NumberTitle', 'off', 'Position', Pos);
	set(gca, 'LineWidth', 1.5, 'FontName', 'Times New Roman', 'FontSize', 20);
	hold on;
	box on;
	axis on;

	%% plot field;
	switch dimPhy
	case 1
		x = sbulkparam.field.spectra;
		y = zeros(length(x), 1);
		c = sbulkparam.field.val;
	case 2
		x = sbulkparam.field.spectra(:,1);
		y = sbulkparam.field.spectra(:,2);
		c = sbulkparam.field.val;
	case 3
		x = sbulkparam.field.spectra(:,1);
		y = sbulkparam.field.spectra(:,2);
		z = sbulkparam.field.spectra(:,3);
		c = sbulkparam.field.val;
	end
	c = c/max(c); % normalize color;
	%% marker size; [0,1];
	szSet = [0.1, 5; 0.2, 10; 0.4, 20; 0.6, 30; 0.8, 40; 1.0, 50];
	sz = c;		szLen = size(szSet);
	for i = 1:1:szLen; sz(sz<=szSet(i,1)) = szSet(i,2); end;
	%% scatter;
	if ( dimPhy == 3 )
		cmap = jet(100); % obtain the color map 'jet';
		for i = 1:1:szLen
			ind = find(c<=szSet(i,1));
			xTmp = x(ind,:);	yTmp = y(ind,:);	zTmp = z(ind,:);
			cTmp = cmap( round(szSet(i,1)*size(cmap,1)), : );
			plot3(xTmp, yTmp, zTmp, 'o', 'MarkerSize', 0.3*szSet(i,2),...
				'MarkerFaceColor', cTmp, 'MarkerEdgeColor', cTmp);
		end
		view(3);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		zlabel('$Z$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(y), max(y), min(z), max(z)]);
	else
		scatter(x, y, 2*sz, c, 'filled', 'Marker', 'o', 'LineWidth', 1.5);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(y), max(y)]);
	end

	%% color;
	colormap('jet');
%	colormap('gray');
	colorbar('FontName', 'Times New Roman', 'FontSize', 20);
	caxis([0, 1]);

	%% saving
	%set(gca, 'LooseInset', [0,0,0,0]);
	figName = sprintf('./fig/%s/bulk%d_field%d.png', sbulkparam.rsltDir,...
				sbulkparam.sflag, plotIter);
	saveas(gca, figName);
end


%% plot the real-space morphology of bulk phases in the interface framework;
function sbulkparam = fn_disp_bulk_proj_density ( sbulkparam )

	%% read data;
	fname = sprintf('./result/%s/bulk%d_FGJP_density.dat',...
				sbulkparam.rsltDir, sbulkparam.sflag);
	if ( exist(fname) == 0 )
		fprintf('\t       %s not exist.\n', fname);
		sbulkparam.projDensity = [];
		return;
	end
	sbulkparam.projDensity = fn_myfile_read ( fname, sbulkparam.dimPhy );

	%% dimension;
	dimPhy = sbulkparam.dimPhy;
	fprintf('\t ====> plot projected density: dimPhy = %d \n', dimPhy);

	%% set configuration of the figure
	figTitle = sprintf('%s:projDensity%d', sbulkparam.rsltDir, sbulkparam.sflag);
	Pos = [150, 250, 1200, 400];
	figure('Name', figTitle, 'NumberTitle', 'off', 'Position', Pos);
	set(gca, 'LineWidth', 2, 'FontName', 'Times New Roman', 'FontSize', 20);
	hold on;

	%% different plotting cases;
	if ( dimPhy == 2 )
		%% spatial grid;
		x = sbulkparam.projDensity.xx;
		y = sbulkparam.projDensity.yy;
		density = sbulkparam.projDensity.val;
		[yy, xx] = meshgrid(y, x);

		%% plotting;
		pcolor(xx, yy, density);
		shading interp;
		box on;
		axis tight;

		%% create colorbar;
		colormap('jet');
%		colormap('gray');
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(y), max(y)]);
	elseif ( dimPhy == 3 )
		%% screen;
		x = sbulkparam.projDensity.xx;
		y = sbulkparam.projDensity.yy;
		z = sbulkparam.projDensity.zz;
		density = sbulkparam.projDensity.val;
		density = real(density);
		minu = min(density(:));
		maxu = max(density(:));
		isoA = minu + 0.7*(maxu-minu);

		%% spatial grid;
		[yy, xx, zz] = meshgrid(y, x, z);
		tol  = 0.2 * [-1, 1, -1, 1, -1, 1];
		axBC = [min(x), max(x), min(y), max(y), min(z), max(z)] + tol;

		%% plot and modify;
		alpA = 0.8;
		patch(isosurface(xx, yy, zz, density, isoA), 'facecolor', [30, 144, 255]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'none');
		patch(isocaps(xx, yy, zz, density, isoA, 'enclose'), 'facecolor', [0, 206, 209]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'interp'); % flat, interp
		daspect([1, 1, 1]);
		camup([1, 0, 0]);
		campos([25, -55, 5]);
		camlight;
		lighting phong;
		box on;
		axis on;
		axis(axBC);
		view([0, 40]);

		%% create colorbar
		colormap('jet');
%		colormap('gray');
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		zlabel('$Z$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
	end

	%% saving
	switch dimPhy
	case 2
		set(gca, 'LooseInset', [0,0,0.02,0.01]);
	case 3
		set(gca, 'LooseInset', [0,0,0,0]);
	end
	figName = sprintf('./fig/%s/bulk%d_FGJP_density.png',...
					sbulkparam.rsltDir, sbulkparam.sflag);
	saveas(gca, figName);
end


%% plot the real-space morphology of bulk phases in the common space;
function sbulkparam = fn_disp_bulk_common_density ( sbulkparam )

	%% read data;
	fname = sprintf('./result/%s/bulk%d_FGJPre_density.dat',...
				sbulkparam.rsltDir, sbulkparam.sflag);
	if ( exist(fname) == 0 )
		fprintf('\t       %s not exist.\n', fname);
		sbulkparam.comDensity = [];
		return;
	end
	sbulkparam.comDensity = fn_myfile_read ( fname, sbulkparam.dimPhy );

	%% dimension;
	dimPhy = sbulkparam.dimPhy;
	fprintf('\t ====> plot density (in the common space): dimPhy = %d \n', dimPhy);

	%% set configuration of the figure
	figTitle = sprintf('%s:comDensity%d', sbulkparam.rsltDir, sbulkparam.sflag);
	Pos = [150, 250, 1200, 400];
	figure('Name', figTitle, 'NumberTitle', 'off', 'Position', Pos);
	set(gca, 'LineWidth', 2, 'FontName', 'Times New Roman', 'FontSize', 20);
	hold on;

	%% different plotting cases;
	if ( dimPhy == 2 )
		%% spatial grid;
		x = sbulkparam.comDensity.xx;
		y = sbulkparam.comDensity.yy;
		density = sbulkparam.comDensity.val;
		[yy, xx] = meshgrid(y, x);

		%% plotting;
		pcolor(xx, yy, density);
		shading interp;
		box on;
		axis tight;

		%% create colorbar;
		colormap('jet');
%		colormap('gray');
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(y), max(y)]);
	elseif ( dimPhy == 3 )
		%% screen;
		x = sbulkparam.comDensity.xx;
		y = sbulkparam.comDensity.yy;
		z = sbulkparam.comDensity.zz;
		density = sbulkparam.comDensity.val;
		density = real(density);
		minu = min(density(:));
		maxu = max(density(:));
		isoA = minu + 0.7*(maxu-minu);

		%% spatial grid;
		[yy, xx, zz] = meshgrid(y, x, z);
		tol  = 0.2 * [-1, 1, -1, 1, -1, 1];
		axBC = [min(x), max(x), min(y), max(y), min(z), max(z)] + tol;

		%% plot and modify;
		alpA = 1.0;
		patch(isosurface(xx, yy, zz, density, isoA), 'facecolor', [30, 144, 255]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'none');
		patch(isocaps(xx, yy, zz, density, isoA, 'enclose'), 'facecolor', [0, 206, 209]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'interp'); % flat, interp
		daspect([1, 1, 1]);
		camup([1, 0, 0]);
		campos([25, -55, 5]);
		camlight;
		lighting phong;
		box on;
		axis on;
		axis(axBC);
		view([0, 40]);

		%% create colorbar
		colormap('jet');
%		colormap('gray');
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		zlabel('$Z$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
	end

	%% saving
	switch dimPhy
	case 2
		set(gca, 'LooseInset', [0,0,0.02,0.01]);
	case 3
		set(gca, 'LooseInset', [0,0,0,0]);
	end
	figName = sprintf('./fig/%s/bulk%d_FGJPre_density.png',...
					sbulkparam.rsltDir, sbulkparam.sflag);
	saveas(gca, figName);
end


%% plot the real-space morphology of connected bulk phases in the common space;
%% based on 'fn_disp_bulk_common_density';
function fn_disp_bulk_connect_density ( sbulkparam1, sbulkparam2, initDist )

	%% dimension;
	dimPhy = sbulkparam1.dimPhy; % == sbulkparam2.dimPhy;
	fprintf('\t ====> plot connected density (in the common space): dimPhy = %d \n', dimPhy);

	%% find the index of 'initDist' to split density matrix of 'sbulkparam1' and 'sbulkparam2';
	x = sbulkparam1.comDensity.xx; % x should be sorted by increasing;
	ind = find ( x >= initDist );
	ind = ind(1);
	density = sbulkparam1.comDensity.val;
	switch dimPhy
	case 2
		density (ind:1:end, : )	   = sbulkparam2.comDensity.val (ind:1:end, : );
	case 3
		density (ind:1:end, :, : ) = sbulkparam2.comDensity.val (ind:1:end, :, : );
	end
	density = real(density);

	%% set configuration of the figure
	figTitle = sprintf('%s:connectDensity', sbulkparam1.rsltDir);
	Pos = [150, 250, 1200, 400];
	figure('Name', figTitle, 'NumberTitle', 'off', 'Position', Pos);
	set(gca, 'LineWidth', 2, 'FontName', 'Times New Roman', 'FontSize', 20);
	hold on;

	%% different plotting cases;
	if ( dimPhy == 2 )
		%% spatial grid;
		y = sbulkparam1.comDensity.yy;
		[yy, xx] = meshgrid(y, x);

		%% plotting;
		pcolor(xx, yy, density);
		shading interp;
		box on;
		axis tight;

		%% create colorbar;
		colormap('jet');
%		colormap('gray');
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(y), max(y)]);
	elseif ( dimPhy == 3 )
		%% screen;
		y = sbulkparam1.comDensity.yy;
		z = sbulkparam1.comDensity.zz;
		minu = min(density(:));
		maxu = max(density(:));
		isoA = minu + 0.7*(maxu-minu);

		%% spatial grid;
		[yy, xx, zz] = meshgrid(y, x, z);
		tol  = 0.2 * [-1, 1, -1, 1, -1, 1];
		axBC = [min(x), max(x), min(y), max(y), min(z), max(z)] + tol;

		%% plot and modify;
		alpA = 1.0;
		patch(isosurface(xx, yy, zz, density, isoA), 'facecolor', [30, 144, 255]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'none');
		patch(isocaps(xx, yy, zz, density, isoA, 'enclose'), 'facecolor', [0, 206, 209]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'interp'); % flat, interp
		daspect([1, 1, 1]);
		camup([1, 0, 0]);
		campos([25, -55, 5]);
		camlight;
		lighting phong;
		box on;
		axis on;
		axis(axBC);
		view([0, 40]);

		%% create colorbar
		colormap('jet');
%		colormap('gray');
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		zlabel('$Z$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
	end

	%% saving
	switch dimPhy
	case 2
		set(gca, 'LooseInset', [0,0,0.02,0.01]);
	case 3
		set(gca, 'LooseInset', [0,0,0,0]);
	end
	figName = sprintf('./fig/%s/bulk_connect_density.png', sbulkparam1.rsltDir);
	saveas(gca, figName);
end


%% obtain bulk parameters;
function ssysparam = fn_obt_system_param ( ssysparam )

	fname = sprintf('./result/%s/parameter_init.dat', ssysparam.rsltDir);

	data = textread(fname, '%s', 'delimiter', '\n');
	for i = 1:1:length(data)
		dataTmp = split(data{i});
		if ( strcmp(dataTmp{1}, 'save_step') )
			ssysparam.save_step = str2num( dataTmp{3} );
		elseif ( strcmp(dataTmp{1}, 'bulk1_dimPhy') )
			ssysparam.dimPhy = str2num( dataTmp{3} );
		elseif ( strcmp(dataTmp{1}, 'x_range') )
			ssysparam.x_range = str2num( dataTmp{3} );
		elseif ( strcmp(dataTmp{1}, 'y_range') )
			ssysparam.y_range = str2num( dataTmp{3} );
		elseif ( strcmp(dataTmp{1}, 'z_range') )
			ssysparam.z_range = str2num( dataTmp{3} );
		elseif ( strcmp(dataTmp{1}, 'GJP_degree') )
			ssysparam.GJP_degree = str2num( dataTmp{3} );
		elseif ( strcmp(dataTmp{1}, 'LGL_num') )
			ssysparam.LGL_num = str2num( dataTmp{3} );
		elseif ( strcmp(dataTmp{1}, 'Four_num') )
			ssysparam.Four_num = str2num( dataTmp{3} );
		end
	end

	fname = sprintf('./result/%s/parameter_opt.dat', ssysparam.rsltDir);
	data = textread(fname, '%s', 'delimiter', '\n');
	for i = 1:1:length(data)
		dataTmp = split(data{i});
		if ( strcmp(dataTmp{1}, 'nprocs') )
			ssysparam.nprocs = str2num( dataTmp{3} );
		elseif ( strcmp(dataTmp{1}, 'com_projmat_size') )
			row = str2num(dataTmp{3});
			col = str2num(dataTmp{4});
		elseif ( strcmp(dataTmp{1}, 'com_projmat_mat') )
			for j = 1:1:row
				buff = split( data{i+j} );
				for k = 1:1:length(buff)
					projMat(j,k) = str2num(buff{k});
				end
			end
		end
	end
	ssysparam.projMat = projMat;
end


%% plot the real-space morphology of bulk phases in the common space;
function ssysparam = fn_disp_system_density ( ssysparam, step, threshold, flag, gbreg )

	%% read data;
	if strcmp(flag, 'read')
		fname = sprintf('./result/%s/sys_density%d.dat', ssysparam.rsltDir, step);
		if ( exist(fname) == 0 )
			fprintf('\t       %s not exist.\n', fname);
			ssysparam.density = [];
			return;
		end
		ssysparam.density = fn_myfile_read ( fname, ssysparam.dimPhy );
	end

	%% dimension;
	dimPhy = ssysparam.dimPhy;
	fprintf('\t ====> plot interface density: step = %d \t dimPhy = %d \n', step, dimPhy);

	%% data;
	x = ssysparam.density.xx;
	y = ssysparam.density.yy;
	density = ssysparam.density.val;
	density = real(density);
	density(density < threshold ) = 0.0;
	minu = min(density(:));
	maxu = max(density(:));
	fprintf('\t       density: [%.15e, %.15e].\n', minu, maxu);
	fprintf('\t       threshold: %.2e\n', threshold);

	%% set configuration of the figure
%	Pos = [150, 250, 1200, 600];
	Pos = [150, 250, 20*ssysparam.x_range, 10*ssysparam.y_range];
	figTitle = sprintf('%s:sysDensity%d', ssysparam.rsltDir, step);
	figure('Name', figTitle, 'NumberTitle', 'off', 'Position', Pos);
	set(gca, 'LineWidth', 2, 'FontName', 'Times New Roman', 'FontSize', 20);
	hold on;

	%% different plotting cases;
	if ( dimPhy == 2 )
		%% spatial grid;
		[yy, xx] = meshgrid(y, x);

		%% plotting;
		pcolor(xx, yy, density);
		shading interp;
		box on;
		axis tight;

		%% create colorbar;
		colormap('jet');
%		colormap('gray');
%		caxis([-4.701846044416446e-01, 1.377654373892305e+00]); % 12fold;
		caxis([-4.646785742801232e-01, 1.229574198931513e+00]); % 10fold;
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(y), max(y)]);

		if ( nargin == 5 )
			x = gbreg(1) * ones(1, length(y));
			plot(x, y, '--', 'LineWidth', 2.0, 'Color', 'w');
			x = gbreg(2) * ones(1, length(y));
			plot(x, y, '--', 'LineWidth', 2.0, 'Color', 'w');
		end

	elseif ( dimPhy == 3 )
		%% spatial grid;
		z = ssysparam.density.zz;
		isoA = minu + 0.7*(maxu-minu);
		[yy, xx, zz] = meshgrid(y, x, z);
		tol  = 1.0 * [-1, 1, -1, 1, -1, 1];
		axBC = [min(x), max(x), min(y), max(y), min(z), max(z)] + tol;

		%% plot and modify;
		alpA = 1.0;
		patch(isosurface(xx, yy, zz, density, isoA), 'facecolor', [30, 144, 255]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'none');
		patch(isocaps(xx, yy, zz, density, isoA, 'enclose'), 'facecolor', [0, 206, 209]/255,...
			'FaceAlpha', alpA, 'edgecolor', 'interp'); % flat, interp
		daspect([1, 1, 1]);
		camup([1, 0, 0]);
		campos([25, -55, 5]);
		camlight;
		lighting phong;
		box on;
		axis on;
		axis(axBC);
		view([0, 40]);

		%% create colorbar
		colormap('jet');
%		colormap('gray');
%		caxis([-4.728842043062894e-01, 1.376037897910467e+00]);
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		zlabel('$Z$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
	end

	%% saving
	switch dimPhy
	case 2
		set(gca, 'LooseInset', [0,0,0.02,0.01]);
	case 3
		set(gca, 'LooseInset', [0,0,0,0]);
	end
	if ( step == -2 )
		figName = sprintf('./fig/%s/recover/recover_density%d.png', ssysparam.rsltDir, step);
	else
		if ( threshold < -100 )
			figName = sprintf('./fig/%s/sys_density%d.png', ssysparam.rsltDir, step);
		else
			thresholdStr = sprintf('%.2f', threshold);	
			for i = 1:1:length(thresholdStr)
				if ( thresholdStr(i) == '.' )
					thresholdStr(i) = 's';
				end
			end
			figName = sprintf('./fig/%s/screen%s_density%d.png', ssysparam.rsltDir, thresholdStr, step);
		end
	end
%	saveas(gca, figName);
	print(gcf, '-dpng', '-r300', figName);
end


%% save screened data;
function ssysparam = fn_save_system_density ( ssysparam, step, threshold )

	%% read data;
	fname = sprintf('./result/%s/sys_density%d.dat', ssysparam.rsltDir, step);
	if ( exist(fname) == 0 )
		fprintf('\t       %s not exist.\n', fname);
		ssysparam.density = [];
		return;
	end
	ssysparam.density = fn_myfile_read ( fname, ssysparam.dimPhy );

	%% dimension;
	dimPhy = ssysparam.dimPhy;
	fprintf('\t ====> save interface density: dimPhy = %d \n', dimPhy);
	fprintf('\t       threshold: %.2e\n', threshold);

	%% file name;
	if ( step == -2 )
		fileName = sprintf('./fig/%s/recover/recover_density%d.txt', ssysparam.rsltDir, step);
	else
		if ( threshold < -100 )
			fileName = sprintf('./fig/%s/sys_density%d.txt', ssysparam.rsltDir, step);
		else
			thresholdStr = sprintf('%.2f', threshold);	
			for i = 1:1:length(thresholdStr)
				if ( thresholdStr(i) == '.' )
					thresholdStr(i) = 's';
				end
			end
			fileName = sprintf('./fig/%s/screen%s_density%d.txt', ssysparam.rsltDir, thresholdStr, step);
		end
	end

	%% data;
	x = ssysparam.density.xx;
	y = ssysparam.density.yy;
	density = ssysparam.density.val;
	density = real(density(:));
	ind = find ( density >= threshold ); % save data;
	len = length(ind);
	fprintf('\t       len: %d\n', len);
	fid = fopen(fileName, 'w');

	%% for ovito
	fprintf(fid, '%d\n', len);
	fprintf(fid, 'Kommentar\n');

	%% different cases;
	if ( dimPhy == 2 )
		%% spatial grid;
		[yy, xx] = meshgrid(y, x);
		xx = xx(:);
		yy = yy(:);
		data = [xx(ind), yy(ind), density(ind)];
		for i = 1:1:len
			fprintf(fid, '%+.6E\t%+.6E\t%+.6E\n', data(i,1), data(i,2), data(i,3));
		end
	elseif ( dimPhy == 3 )
		%% spatial grid;
		z = ssysparam.density.zz;
		[yy, xx, zz] = meshgrid(y, x, z);
		xx = xx(:);
		yy = yy(:);
		zz = zz(:);
		data = [xx(ind), yy(ind), zz(ind), density(ind)];
		for i = 1:1:len
			fprintf(fid, '%+.6E\t%+.6E\t%+.6E\t%+.6E\n', data(i,1), data(i,2), data(i,3), data(i,4));
		end
	end
	fclose(fid);
end




%% save vtk;
function fn_vtk_system_density ( ssysparam, step )

	%% read data;
	if ( isempty(ssysparam.density) )
		fname = sprintf('./result/%s/sys_density%d.dat', ssysparam.rsltDir, step);
		if ( exist(fname) == 0 )
			fprintf('\t       %s not exist.\n', fname);
			return;
		end
		ssysparam.density = fn_myfile_read ( fname, ssysparam.dimPhy );
	end
	if ( step == -2 )
		vtkName = sprintf('./vtk/%s/recover/recover_density%d.vtk', ssysparam.rsltDir, step);
	else
		vtkName = sprintf('./vtk/%s/sys_density%d.vtk', ssysparam.rsltDir, step);
	end

	%% dimension;
	dimPhy = ssysparam.dimPhy;
	fprintf('\t       vtk density (in the interface system): dimPhy = %d \n', dimPhy);

	%% different cases;
	if ( dimPhy == 3 )
		%% screen;
		x = ssysparam.density.xx;
		y = ssysparam.density.yy;
		z = ssysparam.density.zz;
		density = ssysparam.density.val;
		density = real(density);

		%% spatial grid;
		[yy, xx, zz] = meshgrid(y, x, z);
		vtkwrite(vtkName,'STRUCTURED_GRID', xx, yy, zz, 'scalars', 'density', density, 'Precision',13);
	else
		fprintf('\t       vtk files only support 3D.\n');
	end
	fprintf('\n');
end


%% obtain information of interface field (Fourier spectra);
function ssysparam = fn_obt_system_field ( ssysparam, step )

	fprintf('\t ====> obtain interface field: step = %d \n', step);

	%% parameters;
%	nprocs	  = ssysparam.nprocs;
	nprocs	  = 1;
	projMat	  = ssysparam.projMat; % common R'PB;
	[dimRePhy, dimReCpt] = size(projMat);

	%% read 'x';
	fname = sprintf('./result/%s/x.dat', ssysparam.rsltDir);
	fid = fopen(fname, 'r');
	tmp = textscan(fid, '%.15f');
	len = tmp{1}(3);
	x = tmp{1}(4:end);

	%% read 'w';
	fname = sprintf('./result/%s/w.dat', ssysparam.rsltDir);
	fid = fopen(fname, 'r');
	tmp = textscan(fid, '%.15f');
	len = tmp{1}(3);
	w = tmp{1}(4:end);

	%% read 'd0GJP';
	fname = sprintf('./result/%s/d0GJP.dat', ssysparam.rsltDir);
	fid = fopen(fname, 'r');
	tmp = textscan(fid, '%.15f');
	xlen = tmp{1}(1);
	nd = tmp{1}(2);
	len = tmp{1}(3);
	JJ = tmp{1}(4:end);
	JJ = reshape(JJ, [nd, xlen]);

	%% read 'rhoJCplx';
	rhoJCplx = [];
	for myrank = 0:1:nprocs-1
		fname = sprintf('./result/%s/rank%d/sys_rhoJCplx%d.dat',...
				ssysparam.rsltDir, myrank, step);
		if ( exist(fname) == 0 )
			fprintf('\t       %s not exist.\n', fname);
		end
		fid = fopen(fname, 'r');
		tmp = textscan(fid, '%.15f');
		nd  = tmp{1}(1);
		cplxReDofs = tmp{1}(2);
		len = tmp{1}(3);
		val = tmp{1}(4:2:end) + 1i*tmp{1}(5:2:end);
		rhoJCplx = [rhoJCplx; val];
	end
	rhoJCplx = reshape( rhoJCplx, [cplxReDofs, nd]);
	rhoJCplx( isnan(rhoJCplx) ) = 0;

	%% read 'd0bnd';
	d0bnd = [];
	for myrank = 0:1:nprocs-1
		fname = sprintf('./result/%s/rank%d/sys_d0bnd.dat',...
				ssysparam.rsltDir, myrank);
		if ( exist(fname) == 0 )
			fprintf('\t       %s not exist.\n', fname);
		end
		fid = fopen(fname, 'r');
		tmp = textscan(fid, '%.15f');
		xlen  = tmp{1}(1);
		cplxReDofs = tmp{1}(2);
		len = tmp{1}(3);
		val = tmp{1}(4:2:end) + 1i*tmp{1}(5:2:end);
		d0bnd = [d0bnd; val];
	end
	d0bnd = reshape( d0bnd, [cplxReDofs, xlen]);

	%% compute the real-space order parameter by summation
	rhoCplx = rhoJCplx * JJ + d0bnd;

	%% generate the corresponding Fourier spectra;
	Ny = ssysparam.Four_num;
	kk1d = [0:1:Ny/2-1, -Ny/2:1:-1].';
	ncpt = Ny * ones(1, dimReCpt);
	kk = zeros(Ny^dimReCpt, dimReCpt+dimRePhy);
	
	switch dimReCpt
	case 1
		ind = 0;
		for j1 = 1:1:Ny
			ind = ind + 1;
			ktmp = kk1d(j1);
			coordinate = ktmp * projMat.';
			kk(ind,:) = [ktmp, coordinate];
		end
	case 2
		ind = 0;
		for j1 = 1:1:Ny
			for j2 = 1:1:Ny
				ind = ind + 1;
				ktmp = [kk1d(j1), kk1d(j2)];
				coordinate = ktmp * projMat.';
				kk(ind,:) = [ktmp, coordinate];
			end
		end
	case 3
		ind = 0;
		for j1 = 1:1:Ny
			for j2 = 1:1:Ny
				for j3 = 1:1:Ny
					ind = ind + 1;
					ktmp = [kk1d(j1), kk1d(j2), kk1d(j3)];
					coordinate = ktmp * projMat.';
					kk(ind,:) = [ktmp, coordinate];
				end
			end
		end
	case 4
		ind = 0;
		for j1 = 1:1:Ny
			for j2 = 1:1:Ny
				for j3 = 1:1:Ny
					for j4 = 1:1:Ny
						ind = ind + 1;
						ktmp = [kk1d(j1), kk1d(j2), kk1d(j3), kk1d(j4)];
						coordinate = ktmp * projMat.';
						kk(ind,:) = [ktmp, coordinate];
					end
				end
			end
		end
	case 5
		ind = 0;
		for j1 = 1:1:Ny
			for j2 = 1:1:Ny
				for j3 = 1:1:Ny
					for j4 = 1:1:Ny
						for j5 = 1:1:Ny
							ind = ind + 1;
							ktmp = [kk1d(j1), kk1d(j2), kk1d(j3), kk1d(j4), kk1d(j5)];
							coordinate = ktmp * projMat.';
							kk(ind,:) = [ktmp, coordinate];
						end
					end
				end
			end
		end
	case 6
		ind = 0;
		for j1 = 1:1:Ny
			for j2 = 1:1:Ny
				for j3 = 1:1:Ny
					for j4 = 1:1:Ny
						for j5 = 1:1:Ny
							for j6 = 1:1:Ny
								ind = ind + 1;
								ktmp = [kk1d(j1), kk1d(j2), kk1d(j3), kk1d(j4), kk1d(j5), kk1d(j6)];
								coordinate = ktmp * projMat.';
								kk(ind,:) = [ktmp, coordinate];
							end
						end
					end
				end
			end
		end
	end
	
	%% save into the structure body;
	ssysparam.field = struct('kk', kk, 'rhoCplx', rhoCplx);
	ssysparam.base  = struct('x', x, 'w', w,...
		'rhoJCplx', rhoJCplx, 'JJ', JJ, 'd0bnd', d0bnd);
end


%% sort spectra according to the intensity;
function [specIndComm, specIndAll] = fn_sort_system_field ( ssysparam, step, TOL )

	fprintf('\t       sort interface field\n');

	%% parameters;
	projMat	  = ssysparam.projMat; % common R'PB;
	[dimRePhy, dimReCpt] = size(projMat);
	kk = ssysparam.field.kk;
	rhoCplx = ssysparam.field.rhoCplx;
	[nr, nc] = size(rhoCplx);

	%% generate spectra;
	specIndComm = [];
	filename = sprintf('./fig/%s/sys_field%d.dat', ssysparam.rsltDir, step);
	fid = fopen(filename, 'w'); 
	for xslice = 1:1:nc
		fprintf(fid, 'xslice = %d\n', xslice);
		[val, ind] = sort( abs(rhoCplx(:,xslice)), 'descend');
		xspecInd = find( val>TOL );
		%% screen spectra indicates;
		if ( xslice == 1 )
			specIndComm = ind(xspecInd); % save common indicates;
			specIndAll = ind(xspecInd); % save all appeared indicates;
		else
			specIndComm = intersect(specIndComm, ind(xspecInd));
			specIndAll = union(specIndAll, ind(xspecInd));
		end
		len = length( xspecInd );
		for i = 1:1:len
			fprintf(fid, 'ind%*d:\t', 6, ind(i));
			for j = 1:1:dimReCpt
				fprintf(fid, ' % d', kk(ind(i),j));
			end
			fprintf(fid, '\t');
			for j = dimReCpt+1:1:dimReCpt+dimRePhy
				fprintf(fid, ' % .6f', kk(ind(i),j));
			end
			tmp = rhoCplx(ind(i),xslice);
			fprintf(fid, '\t --> %.6f,%.6f\n', real(tmp), imag(tmp));
		end
		fprintf(fid, '\n');
	end
	fprintf(fid, '\n\n');
	fprintf(fid, 'specIndComm: \n');
	for i = 1:1:length(specIndComm)
		fprintf(fid, 'ind%*d:\t', 6, specIndComm(i));
		for j = 1:1:dimReCpt
			fprintf(fid, ' % d', kk(specIndComm(i),j));
		end
		fprintf(fid, '\t');
		for j = dimReCpt+1:1:dimReCpt+dimRePhy
			fprintf(fid, ' % .6f', kk(specIndComm(i),j));
		end
		fprintf(fid, '\n');
	end
	fprintf(fid, '\n');
	fprintf(fid, 'specIndAll: \n');
	for i = 1:1:length(specIndAll)
		fprintf(fid, 'ind%*d:\t', 6, specIndAll(i));
		for j = 1:1:dimReCpt
			fprintf(fid, ' % d', kk(specIndAll(i),j));
		end
		fprintf(fid, '\t');
		for j = dimReCpt+1:1:dimReCpt+dimRePhy
			fprintf(fid, ' % .6f', kk(specIndAll(i),j));
		end
		fprintf(fid, '\n');
	end
	fprintf(fid, '\n\n');
	fclose(fid);

	fprintf('specIndComm: ');
	for i = 1:1:length(specIndComm)
		fprintf(' %d', specIndComm(i));
	end
	fprintf('\n');
	fprintf('specIndAll : ');
	for i = 1:1:length(specIndAll)
		fprintf(' %d', specIndAll(i));
	end
	fprintf('\n');
end


%% plot the slices of interface field (Fourier spectra);
function fn_disp_system_field_xslice ( ssysparam, step, xslice, TOL )

	%% parameters;
%	nprocs	  = ssysparam.nprocs;
	nprocs	  = 1;
	projMat	  = ssysparam.projMat; % common R'PB;
	[dimRePhy, dimReCpt] = size(projMat);
	kk = ssysparam.field.kk;
	rhoCplx = ssysparam.field.rhoCplx;

	%% generate spectra;
	spectra = [];
	FourCoeff = [];
	[nr, nc] = size(kk);
	for i = 1:1:nr
		val = abs(rhoCplx(i, xslice));
		if ( val > TOL )
			spectra(end+1,:) = kk(i,:);
			FourCoeff(end+1) = val;
		end
	end
	[val, ind]  = sort(FourCoeff, 'descend');
	spectra = spectra(ind,:);

	fprintf('\t ====> plot interface slices:  ');
	fprintf('dimReCpt = %d, step = %d, xslice = %d \n', dimReCpt, step, xslice);

	%% set configuration of the figure
	Pos = [200, 250, 800, 600];
	figTitle = sprintf('%s:field%d_xslice%d', ssysparam.rsltDir, step, xslice);
	figure('Name', figTitle, 'NumberTitle', 'off', 'Position', Pos);
	set(gca, 'LineWidth', 1.5, 'FontName', 'Times New Roman', 'FontSize', 20);
	hold on;
	box on;
	axis on;

	%% plot field;
	switch dimRePhy
	case 1
		x = spectra(:,dimReCpt+1);
	case 2
		x = spectra(:,dimReCpt+1);
		y = spectra(:,dimReCpt+2);
	case 3
		x = spectra(:,dimReCpt+1);
		y = spectra(:,dimReCpt+2);
		z = spectra(:,dimReCpt+3);
	end
	c = val/max(abs(val)); % normalize color;
	%% marker size; [0,1];
	szSet = [0.1, 5; 0.2, 10; 0.4, 20; 0.6, 30; 0.8, 40; 1.0, 50];
	sz = c;		szLen = size(szSet);
	for i = 1:1:szLen; sz(sz<=szSet(i,1)) = szSet(i,2); end;
	%% scatter;
	switch dimRePhy
	case 1
		plot(x, val, 'o', 'MarkerSize', 10,...
		   	'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm');
		xlabel('$\tilde{P}\mathbf{h}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$|\hat{\phi}(x,\mathbf{h})|$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(val), max(val)]);
	case 2
		scatter(x, y, 2*sz, c, 'filled', 'Marker', 'o', 'LineWidth', 1.5);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(y), max(y)]);
	case 3
		cmap = jet(100); % obtain the color map 'jet';
		for i = 1:1:szLen
			ind = find(c<=szSet(i,1));
			xTmp = x(ind,:);	yTmp = y(ind,:);	zTmp = z(ind,:);
			cTmp = cmap( round(szSet(i,1)*size(cmap,1)), : );
			plot3(xTmp, yTmp, zTmp, 'o', 'MarkerSize', 0.3*szSet(i,2),...
				'MarkerFaceColor', cTmp, 'MarkerEdgeColor', cTmp);
		end
		view(3);
		xlabel('$X$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		ylabel('$Y$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		zlabel('$Z$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
		axis([min(x), max(x), min(y), max(y), min(z), max(z)]);

		%% color;
		colormap('jet');
	%	colormap('gray');
		colorbar('FontName', 'Times New Roman', 'FontSize', 20);
		caxis([0, 1]);
	end

	%% saving
	%set(gca, 'LooseInset', [0,0,0,0]);
	figName = sprintf('./fig/%s/xslice/sys_field%d_xslice%d.png',...
			ssysparam.rsltDir, step, xslice);
	saveas(gca, figName);
end


%% plot the curves of interface field against x;
function fn_disp_system_field_xcurve ( ssysparam, step, specInd )

	%% parameters;
%	nprocs	  = ssysparam.nprocs;
	nprocs	  = 1;
	x_range	  = ssysparam.x_range;
	y_range	  = ssysparam.y_range;
	projMat	  = ssysparam.projMat; % common R'PB;
	[dimRePhy, dimReCpt] = size(projMat);

	%% field;
	kk = ssysparam.field.kk;
	rhoCplx = ssysparam.field.rhoCplx;
	y = abs( rhoCplx(specInd, :) );
	x = linspace(-x_range, x_range, length(y));

	fprintf('\t ====> plot interface curves:  ');
	fprintf('x_range = %.2f, step = %d, specInd = %d \n', x_range, step, specInd);

	%% set configuration of the figure
	Pos = [200, 250, 1000, 400];
	figTitle = sprintf('%s:field%d_specInd%d', ssysparam.rsltDir, step, specInd);
	figure('Name', figTitle, 'NumberTitle', 'off', 'Position', Pos);
	set(gca, 'LineWidth', 1.5, 'FontName', 'Times New Roman', 'FontSize', 20);
	hold on;
	box on;
	axis on;

	%% plot field;
	plot(x, y, '-o', 'Color', 'b', 'LineWidth', 1.5, 'MarkerSize', 5,...
		'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'm');
	xlabel('$x$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
	ylabel('$|\hat{\phi}(x,\mathbf{h})|$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 20);
	axis([min(x), max(x), min(y), max(y)]);

	%% saving
	set(gca, 'LooseInset', [0,0,0.01,0.02]);
	figName = sprintf('./fig/%s/xcurve/sys_field%d_specInd%d.png',...
			ssysparam.rsltDir, step, specInd);
	saveas(gca, figName);
end


%% recover interface density profile by the screened spectra;
function fn_obt_system_recover_density ( ssysparam, step, TOL )

	fprintf('\t       recover rhoJCplx\n');

	%% load data;
	rhoCplx	 = ssysparam.field.rhoCplx;
	rhoJCplx = ssysparam.base.rhoJCplx;
	JJ		 = ssysparam.base.JJ;
	d0bnd	 = ssysparam.base.d0bnd;
	x		 = ssysparam.base.x;
	w		 = ssysparam.base.w;
	[nr, nc] = size(rhoJCplx);
	len = nr * nc;

	%% screen;
%	rhoCplx ( abs(rhoCplx) <= TOL ) = 0.0;

	%% recover;
	mat0 = JJ * diag(w) * JJ.';
	mat0 = kron(speye(nr), mat0);
	mat0 = sparse(mat0.*(abs(mat0)>10^(-8)));
	rhs = JJ * diag(w) * ( rhoCplx - d0bnd ).';
	rhoJCplx = mat0 \ rhs(:);
	rhoJCplx( isnan(rhoJCplx) ) = 0;

	%% save;
	rhoJCplx = reshape( rhoJCplx, [nc, nr] ).';
	rhoJCplx = rhoJCplx(:);
	fname = sprintf('./fig/%s/recover%d.dat', ssysparam.rsltDir, step);
	fid = fopen(fname, 'w');
	fprintf(fid, '%d\t%d\t%d\n', nc, nr, len);
	for i = 1:1:len
		fprintf(fid, "%+.15E\t%+.15E\n", real(rhoJCplx(i)), imag(rhoJCplx(i)));
	end
	fclose(fid);
end


%% file read;
function srslt = fn_myfile_read ( fname, dim )

	fid = fopen(fname, 'r');
	tmp = textscan(fid, '%.15f');

	%% obtain grids and density;
	switch dim
	case 1
		nx = tmp{1}(1);		% the number of discrete x;
		nn = tmp{1}(2);		% the number of total data;
		xx = tmp{1}(3:1:2+nx);
		val = tmp{1}(3+nx:1:2+nx+nn);
		srslt = struct('xx', xx, 'val', val);
	case 2
		nx = tmp{1}(1);		% the number of discrete x;
		ny = tmp{1}(2);		% the number of discrete y;
		nn = tmp{1}(3);		% the number of total data;
		xx = tmp{1}(4:1:3+nx);
		yy = tmp{1}(4+nx:1:3+nx+ny);
		val = tmp{1}(4+nx+ny:1:3+nx+ny+nn);
		val = permute( reshape(val, [ny, nx]), [2,1] );
		srslt = struct('xx', xx, 'yy', yy, 'val', val);
	case 3
		nx = tmp{1}(1);		% the number of discrete x;
		ny = tmp{1}(2);		% the number of discrete y;
		nz = tmp{1}(3);		% the number of discrete z;
		nn = tmp{1}(4);		% the number of total data;
		xx = tmp{1}(5:1:4+nx);
		yy = tmp{1}(5+nx:1:4+nx+ny);
		zz = tmp{1}(5+nx+ny:1:4+nx+ny+nz);
		val = tmp{1}(5+nx+ny+nz:1:4+nx+ny+nz+nn);
		val = permute( reshape(val, [nz, ny, nx]), [3,2,1] );
		srslt = struct('xx', xx, 'yy', yy, 'zz', zz, 'val', val);
	end
	fclose(fid);
end


%% matrix expand;
function uc = mymap ( uc_cs, uc_rf, mapflag )
	
	%% matrix size;
	n_cs = size(uc_cs);
	n_rf = size(uc_rf);
	kk = n_rf/2 - n_cs/2;

	%% shift;
	tmp_rf = fftshift(uc_rf);
	tmp_cs = fftshift(uc_cs);

	%% map;
	if	strcmp('cs2rf', mapflag) % c2r;
		switch length(n_cs)
		case 1
			tmp_rf(kk(1)+1:1:kk(1)+n_cs(1)) = tmp_cs; 
		case 2
			tmp_rf(kk(1)+1:1:kk(1)+n_cs(1), kk(2)+1:1:kk(2)+n_cs(2)) = tmp_cs; 
		case 3
			tmp_rf(kk(1)+1:1:kk(1)+n_cs(1), kk(2)+1:1:kk(2)+n_cs(2),...
				kk(3)+1:1:kk(3)+n_cs(3)) = tmp_cs; 
		case 4
			tmp_rf(kk(1)+1:1:kk(1)+n_cs(1), kk(2)+1:1:kk(2)+n_cs(2),...
				kk(3)+1:1:kk(3)+n_cs(3), kk(4)+1:1:kk(4)+n_cs(4)) = tmp_cs; 
		end
		uc = fftshift(tmp_rf);
	elseif strcmp('rf2cs', mapflag) % r2c;
		switch length(n_cs)
		case 1
			tmp_cs = tmp_rf(kk(1)+1:1:kk(1)+n_cs(1)); 
		case 2
			tmp_cs = tmp_rf(kk(1)+1:1:kk(1)+n_cs(1), kk(2)+1:1:kk(2)+n_cs(2)); 
		case 3
			tmp_cs = tmp_rf(kk(1)+1:1:kk(1)+n_cs(1), kk(2)+1:1:kk(2)+n_cs(2),...
				kk(3)+1:1:kk(3)+n_cs(3)); 
		case 4
			tmp_cs = tmp_rf(kk(1)+1:1:kk(1)+n_cs(1), kk(2)+1:1:kk(2)+n_cs(2),...
				kk(3)+1:1:kk(3)+n_cs(3), kk(4)+1:1:kk(4)+n_cs(4)); 
		end
		uc = fftshift(tmp_cs);
	end
end


%% \brief	 draw lattice of bulk phase;
%% \param	 sbulkparam:	the structure body provides parameters;
%% \param	 axBC:			the boundaries of axes; (should be large enough);
function fn_disp_bulk_lattice ( sbulkparam, axBC, plotIter )

	dirBox = 2*pi*inv(sbulkparam.rcpBox);
	rotateMat = sbulkparam.rotateMat';
	if ( plotIter == -2 || plotIter == -3 )
		translVec = sbulkparam.translVec';
	else
		translVec = zeros(1,3);
	end
	if ( plotIter == -3 )
		rotateMat = eye(3);
	end

	%% generate lattices;
	if isdiag(dirBox)
		xblen = dirBox(1,1);
		yblen = dirBox(2,2);
		zblen = dirBox(3,3);
		x = [0, -xblen:-xblen:axBC(1), xblen:xblen:axBC(2)];
		y = [0, -yblen:-yblen:axBC(3), yblen:yblen:axBC(4)];
		z = [0, -zblen:-zblen:axBC(5), zblen:zblen:axBC(6)];
		xlen = length(x);
		ylen = length(y);
		zlen = length(z);
		lattice = [];
		for i = 1:1:xlen
			for j = 1:1:ylen
				coordSta = [x(i), y(j), min(z)] * rotateMat.';
				coordEnd = [x(i), y(j), min(z)] * rotateMat.';
				lattice(end+1,:) = [coordSta, coordEnd];
			end
		end
		for i = 1:1:xlen
			for j = 1:1:zlen
				coordSta = [x(i), min(y), z(j)] * rotateMat.';
				coordEnd = [x(i), max(y), z(j)] * rotateMat.';
				lattice(end+1,:) = [coordSta, coordEnd];
			end
		end
		for i = 1:1:ylen
			for j = 1:1:zlen
				coordSta = [min(x), y(i), z(j)] * rotateMat.';
				coordEnd = [max(x), y(i), z(j)] * rotateMat.';
				lattice(end+1,:) = [coordSta, coordEnd];
			end
		end
	end
	lattice = lattice + [translVec, translVec];
	line(lattice(:,[1,4])', lattice(:,[2,5])', lattice(:,[3,6])',...
		'LineWidth', 2.0, 'Color', 'r');
end


function vtkwrite( filename,dataType,varargin )

	% VTKWRITE Writes 3D Matlab array into VTK file format.
	%  vtkwrite(filename,'structured_grid',x,y,z,'vectors',title,u,v,w) writes
	%  a structured 3D vector data into VTK file, with name specified by the string
	%  filename. (u,v,w) are the vector components at the points (x,y,z). x,y,z
	%  should be 3-D matrices like those generated by meshgrid, where
	%  point(ijk) is specified by x(i,j,k), y(i,j,k) and z(i,j,k).
	%  The matrices x,y,z,u,v,w must all be the same size and contain
	%  corrresponding position and vector component. The string title specifies
	%  the name of the vector field to be saved. 
	%
	%  vtkwrite(filename,'structured_grid',x,y,z,'scalars',title,r) writes a 3D
	%  scalar data into VTK file whose name is specified by the string
	%  filename. r is the scalar value at the points (x,y,z). The matrices
	%  x,y,z,r must all be the same size and contain the corresponding position
	%  and scalar values. 
	%
	%  vtkwrite(filename,'structured_grid',x,y,z,'vectors',title,u,v,w,'scalars',
	%  title2,r) writes a 3D structured grid that contains both vector and scalar values.
	%  x,y,z,u,v,w,r must all be the same size and contain the corresponding
	%  positon, vector and scalar values.
	%
	%  vtkwrite(filename, 'structured_points', title, m) saves matrix m (could
	%  be 1D, 2D or 3D array) into vtk as structured points.
	%
	%  vtkwrite(filename, 'structured_points', title, m, 'spacing', sx, sy, sz)
	%  allows user to specify spacing. (default: 1, 1, 1). This is the aspect
	%  ratio of a single voxel. 
	%
	%  vtkwrite(filename, 'structured_points', title, m, 'origin', ox, oy, oz)
	%  allows user to speicify origin of dataset. (default: 0, 0, 0).
	%
	%  vtkwrite(filename,'unstructured_grid',x,y,z,'vectors',title,u,v,w,'scalars',
	%  title2,r) writes a 3D unstructured grid that contains both vector and scalar values.
	%  x,y,z,u,v,w,r must all be the same size and contain the corresponding
	%  positon, vector and scalar values.
	%  
	%  vtkwrite(filename, 'polydata', 'lines', x, y, z) exports a 3D line where
	%  x,y,z are coordinates of the points that make the line. x, y, z are
	%  vectors containing the coordinates of points of the line, where point(n)
	%  is specified by x(n), y(n) and z(n).
	%
	%  vtkwrite(filename,'polydata','lines',x,y,z,'Precision',n) allows you to
	%  specify precision of the exported number up to n digits after decimal
	%  point. Default precision is 3 digits. 
	%
	%  vtkwrite(filename,'polydata','triangle',x,y,z,tri) exports a list of
	%  triangles where x,y,z are the coordinates of the points and tri is an
	%  m*3 matrix whose rows denote the points of the individual triangles.
	%
	%  vtkwrite(filename,'polydata','tetrahedron',x,y,z,tetra) exports a list
	%  of tetrahedrons where x,y,z are the coordinates of the points
	%  and tetra is an m*4 matrix whose rows denote the points of individual
	%  tetrahedrons. 
	%  
	%  vtkwrite('execute','polydata','lines',x,y,z) will save data with default
	%  filename ''matlab_export.vtk' and automatically loads data into
	%  ParaView. 
	%  
	%  Version 2.3
	%  Copyright, Chaoyuan Yeh, 2016
	%  Codes are modified from William Thielicke and David Gingras's submission.    

	if strcmpi(filename,'execute'), filename = 'matlab_export.vtk'; end
	fid = fopen(filename, 'w'); 
	% VTK files contain five major parts
	% 1. VTK DataFile Version
	fprintf(fid, '# vtk DataFile Version 2.0\n');
	% 2. Title
	fprintf(fid, 'VTK from Matlab\n');


	binaryflag = any(strcmpi(varargin, 'BINARY'));
	if any(strcmpi(varargin, 'PRECISION'))
		precision = num2str(varargin{find(strcmpi(varargin, 'PRECISION'))+1});
	else
		precision = '2';
	end

	switch upper(dataType)
		case 'STRUCTURED_POINTS'
			title = varargin{1};
			m = varargin{2};
			if any(strcmpi(varargin, 'spacing'))
				sx = varargin{find(strcmpi(varargin, 'spacing'))+1};
				sy = varargin{find(strcmpi(varargin, 'spacing'))+2};
				sz = varargin{find(strcmpi(varargin, 'spacing'))+3};
			else
				sx = 1;
				sy = 1;
				sz = 1;
			end
			if any(strcmpi(varargin, 'origin'))
				ox = varargin{find(strcmpi(varargin, 'origin'))+1};
				oy = varargin{find(strcmpi(varargin, 'origin'))+2};
				oz = varargin{find(strcmpi(varargin, 'origin'))+3};
			else
				ox = 0;
				oy = 0;
				oz = 0;
			end
			[nx, ny, nz] = size(m);
			setdataformat(fid, binaryflag);
			
			fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
			fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
			fprintf(fid, ['SPACING ', num2str(sx), ' ', num2str(sy), ' ',...
				num2str(sz), '\n']);
			fprintf(fid, ['ORIGIN ', num2str(ox), ' ', num2str(oy), ' ',...
				num2str(oz), '\n']); 
			fprintf(fid, 'POINT_DATA %d\n', nx*ny*nz);
			fprintf(fid, ['SCALARS ', title, ' float 1\n']);
			fprintf(fid,'LOOKUP_TABLE default\n');
			if ~binaryflag 
				spec = ['%0.', precision, 'f '];
				fprintf(fid, spec, m(:)');
			else
				fwrite(fid, m(:)', 'float', 'b');
			end
			
		case {'STRUCTURED_GRID','UNSTRUCTURED_GRID'}
			% 3. The format data proper is saved in (ASCII or Binary). Use
			% fprintf to write data in the case of ASCII and fwrite for binary.
			if numel(varargin)<6, error('Not enough input arguments'); end
			setdataformat(fid, binaryflag);
	%         fprintf(fid, 'BINARY\n');
			x = varargin{1};
			y = varargin{2};
			z = varargin{3};
			if sum(size(x)==size(y) & size(y)==size(z))~=length(size(x))
				error('Input dimesions do not match')
			end
			n_elements = numel(x);
			% 4. Type of Dataset ( can be STRUCTURED_POINTS, STRUCTURED_GRID,
			% UNSTRUCTURED_GRID, POLYDATA, RECTILINEAR_GRID or FIELD )
			% This part, dataset structure, begins with a line containing the
			% keyword 'DATASET' followed by a keyword describing the type of dataset.
			% Then the geomettry part describes geometry and topology of the dataset.
			if strcmpi(dataType,'STRUCTURED_GRID')
				fprintf(fid, 'DATASET STRUCTURED_GRID\n');
				fprintf(fid, 'DIMENSIONS %d %d %d\n', size(x,1), size(x,2), size(x,3));
			else
				fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
			end
			fprintf(fid, ['POINTS ' num2str(n_elements) ' float\n']);
			output = [x(:)'; y(:)'; z(:)'];
			
			if ~binaryflag
				spec = ['%0.', precision, 'f '];
				fprintf(fid, spec, output);
			else
				fwrite(fid, output, 'float', 'b');
			end
			% 5.This final part describe the dataset attributes and begins with the
			% keywords 'POINT_DATA' or 'CELL_DATA', followed by an integer number
			% specifying the number of points of cells. Other keyword/data combination
			% then define the actual dataset attribute values.
			fprintf(fid, ['\nPOINT_DATA ' num2str(n_elements)]);
			% Parse remaining argument.
			vidx = find(strcmpi(varargin,'VECTORS'));
			sidx = find(strcmpi(varargin,'SCALARS'));
			if vidx~=0
				for ii = 1:length(vidx)
					title = varargin{vidx(ii)+1};
					% Data enteries begin with a keyword specifying data type
					% and numeric format.
					fprintf(fid, ['\nVECTORS ', title,' float\n']);
					output = [varargin{ vidx(ii) + 2 }(:)';...
							  varargin{ vidx(ii) + 3 }(:)';...
							  varargin{ vidx(ii) + 4 }(:)'];

					if ~binaryflag
						spec = ['%0.', precision, 'f '];
						fprintf(fid, spec, output);
					else
						fwrite(fid, output, 'float', 'b');
					end
	%                 fwrite(fid, [reshape(varargin{vidx(ii)+2},1,n_elements);...
	%                 reshape(varargin{vidx(ii)+3},1,n_elements);...
	%                 reshape(varargin{vidx(ii)+4},1,n_elements)],'float','b');
				end
			end
			if sidx~=0
				for ii = 1:length(sidx)
					title = varargin{sidx(ii)+1};
					fprintf(fid, ['\nSCALARS ', title,' float\n']);
					fprintf(fid, 'LOOKUP_TABLE default\n');
					if ~binaryflag
						spec = ['%0.', precision, 'f '];
						fprintf(fid, spec, varargin{ sidx(ii) + 2});
					else
						fwrite(fid, varargin{ sidx(ii) + 2}, 'float', 'b');
					end
	%                 fwrite(fid, reshape(varargin{sidx(ii)+2},1,n_elements),'float','b');
				end
			end
			
		case 'POLYDATA'

			fprintf(fid, 'ASCII\n');
			if numel(varargin)<4, error('Not enough input arguments'); end
			x = varargin{2}(:);
			y = varargin{3}(:);
			z = varargin{4}(:);
			if numel(varargin)<4, error('Not enough input arguments'); end
			if sum(size(x)==size(y) & size(y)==size(z))~= length(size(x))
				error('Input dimesions do not match')
			end
			n_elements = numel(x);
			fprintf(fid, 'DATASET POLYDATA\n');
			if mod(n_elements,3)==1
				x(n_elements+1:n_elements+2,1)=[0;0];
				y(n_elements+1:n_elements+2,1)=[0;0];
				z(n_elements+1:n_elements+2,1)=[0;0];
			elseif mod(n_elements,3)==2
				x(n_elements+1,1)=0;
				y(n_elements+1,1)=0;
				z(n_elements+1,1)=0;
			end
			nbpoint = numel(x);
			fprintf(fid, ['POINTS ' num2str(nbpoint) ' float\n']);
			
			spec = [repmat(['%0.', precision, 'f '], 1, 9), '\n'];
			
			output = [x(1:3:end-2), y(1:3:end-2), z(1:3:end-2),...
					  x(2:3:end-1), y(2:3:end-1), z(2:3:end-1),...
					  x(3:3:end), y(3:3:end), z(3:3:end)]';
				  
			fprintf(fid, spec, output);
			
			switch upper(varargin{1})
				case 'LINES'
					if mod(n_elements,2)==0
						nbLine = 2*n_elements-2;
					else
						nbLine = 2*(n_elements-1);
					end
					conn1 = zeros(nbLine,1);
					conn2 = zeros(nbLine,1);
					conn2(1:nbLine/2) = 1:nbLine/2;
					conn1(1:nbLine/2) = conn2(1:nbLine/2)-1;
					conn1(nbLine/2+1:end) = 1:nbLine/2;
					conn2(nbLine/2+1:end) = conn1(nbLine/2+1:end)-1;
					fprintf(fid,'\nLINES %d %d\n',nbLine,3*nbLine);
					fprintf(fid,'2 %d %d\n',[conn1';conn2']);
				case 'TRIANGLE'
					ntri = length(varargin{5});
					fprintf(fid,'\nPOLYGONS %d %d\n',ntri,4*ntri);
					fprintf(fid,'3 %d %d %d\n',(varargin{5}-1)');
				case 'TETRAHEDRON'
					ntetra = length(varargin{5});
					fprintf(fid,'\nPOLYGONS %d %d\n',ntetra,5*ntetra);
					fprintf(fid,'4 %d %d %d %d\n',(varargin{5}-1)');
			end     
	end
	fclose(fid);
	if strcmpi(filename,'matlab_export.vtk')
		switch computer
			case {'PCWIN','PCWIN64'}
				!paraview.exe --data='matlab_export.vtk' &
				% Exclamation point character is a shell escape, the rest of the
				% input line will be sent to operating system. It can not take
				% variables, though. The & at the end of line will return control to 
				% Matlab even when the outside process is still running. 
			case {'GLNXA64','MACI64'}
				!paraview --data='matlab_export.vtk' &
		end
	end
	end

	function setdataformat(fid, binaryflag)

	if ~binaryflag
		fprintf(fid, 'ASCII\n');
	else
		fprintf(fid, 'BINARY\n');
	end
end

