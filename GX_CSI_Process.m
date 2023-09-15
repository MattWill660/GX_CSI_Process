function status = GX_CSI_Process(GX_file, H_file)

    path = fileparts(GX_file);

    [data.gx, info.gx] = load_gx(GX_file); % t, kx, ky, kz - load gx data
    info.gx.gen = get_gen_list_info(GX_file);

    [data.h, info.h] = load_H(H_file); % kx, ky, kz, c - load h data
    info.h.gen = get_gen_list_info(H_file);
    
    info = get_acq_info(info); % load necessary acquisition info for processing
    
    [img_gx_fids] = spatial_fft(data.gx); % transform to t, x, y, z

    [fit_gx_mean_fid] = fit_fid(mean(img_gx_fids, [2 3 4]), info.gx, []); % fit mean to get initial guesses

    img_gx_fits(info.gx.nks(1), info.gx.nks(2), info.gx.nks(3)) = {fit_gx_mean_fid}; % create cell aray to hold voxel fits

    % fit all spectra
    for x = 1:info.gx.nks(1)
        for y = 1:info.gx.nks(2)
            parfor z = 1:info.gx.nks(3)
                fit = fit_fid(img_gx_fids(:,x,y,z), info.gx, fit_gx_mean_fid.new_guesses);
                img_gx_fits(x,y,z) = {fit};
                disp(['x/y/z = ', num2str(x), '/', num2str(y), '/', num2str(z),' of ', num2str(info.gx.nks(1)),'/',num2str(info.gx.nks(2)),'/',num2str(info.gx.nks(3))])
            end
        end
    end

    % reconstruct H
    images.h = h_recon(data.h, info.h);

    % reconstruct xe
    images.xe = decomp_fits(img_gx_fits);

    % split / calculate quantitative maps
    [maps, images] = calc_maps(images);

    % interp to match
    [images, maps] = interp_images(images, maps);

    % calculate quantitative results
    [results, maps] = quantify_results(maps, images);

    % create figures
    hlims = [prctile(abs(images.h(:)),3) prctile(abs(images.h(:)),97)];
    % qualitative
    figures.vent = make_fig(abs(images.xe.gas),'hot',...
        [],results.amap.gas,...
        abs(images.h),[],hlims);
    figures.mem = make_fig(abs(images.xe.membrane),'hot',...
        [],results.amap.mem,...
        abs(images.h),[],hlims);
    figures.rbc = make_fig(abs(images.xe.rbc),'hot',...
        [],results.amap.rbc,...
        abs(images.h),[],hlims);
    % quant uptake
    figures.mem_uptake = make_fig(maps.mem_quant,'hot',...
        [0 0.038],results.amap.mem,...
        abs(images.h),[],hlims);
    figures.rbc_transfer = make_fig(maps.rbc_quant,'hot',...
        [0 0.012],results.amap.rbc,...
        abs(images.h),[],hlims);
    figures.rbc2mem = make_fig(maps.rbc2mem_quant,'hot',...
        [0 1],results.amap.mem,...
        abs(images.h),[],hlims);
    % quant frequency
    figures.gas_freq = make_fig(maps.gas_f,'hot',...
        [-7400 -7300],results.amap.gas,...
        abs(images.h),[],hlims);
    figures.membrane_freq = make_fig(maps.membrane_f-results.stats.gas_f.median,'hot',...
        [6900 7000],results.amap.mem,...
        abs(images.h),[],hlims);
    figures.rbc_freq = make_fig(maps.rbc_f-results.stats.gas_f.median,'hot',...
        [7650 7750],results.amap.rbc,...
        abs(images.h),[],hlims);
    % quant t2
    figures.gas_t2 = make_fig(maps.gas_t2,'hot',...
        [0 20],results.amap.gas,...
        abs(images.h),[],hlims);
    figures.membrane_t2 = make_fig(maps.membrane_t2,'hot',...
        [0 2],results.amap.mem,...
        abs(images.h),[],hlims);
    figures.rbc_t2 = make_fig(maps.rbc_t2,'hot',...
        [0 2],results.amap.rbc,...
        abs(images.h),[],hlims);

    % create report  % not implemented atm; only worry about if used more often
%     rep = make_report(images, maps, figures, results);

    % save results
    save_images(images,path);
    save_maps(maps,path);
    save_figures(figures,path);
%     save_report(rep,path); % not implemented atm; only worry about if used more often

    % report results
    status = 0;

    % save mat 
    save(fullfile(path, 'gx_csi_process.mat'));
end


function [data, info] = load_gx(GX_file)
    % only data list option currently
    [data, list_info] = loadLISTDATAspec(GX_file);
    data = squeeze(data);
    info.list_info = list_info;
    info.nt = size(data,1);
    info.nks = [size(data,2), size(data,3), size(data,4)];
    info.nc = 1;
end

function [data, info] = load_H(H_file)
    % only data list option currently
    [data, list_info] = loadLISTDATA(H_file);
    data = squeeze(data);
    info.list_info = list_info;
    info.nt = 1;
    info.nks = [size(data,1), size(data,2), size(data,3)];
    info.nc = size(data,4);
end

function img_fids = spatial_fft(gx_data)
    img_fids = fftshift(fft(fftshift(gx_data,2),[],2),2); % x
    img_fids = fftshift(fft(fftshift(img_fids,3),[],3),3); % y
    img_fids = fftshift(fft(fftshift(img_fids,4),[],4),4); % z
end

function info = get_acq_info(info)
    % how to automate?
    info.gx.dwell = 62.5 /1000000; % us to s
    info.h.dwell = 62.5 /1000000; % us to s % doesnt matter atm
end

function fit_result = fit_fid(fid, info, guesses)
    
    if isempty(info)
        info.dwell = 62.5 /1000000; %us to s
        info.nt = size(fid,1);
    end
    
    if isempty(guesses)
        % initial guesses rbc/mem/gas
        guesses.area_guess = [abs(fid(1)) abs(fid(1)) abs(fid(1))]; 
        guesses.area_lowerBounds = [0 0 0]; % Always positive
        guesses.area_upperBounds = 2*guesses.area_guess;
        guesses.freq_guess = [353 -353 -7350]; % assume 3T and receive between dissolved peaks
        guesses.freq_lowerBounds = [0 -2000 -9350];
        guesses.freq_upperBounds = [2000 0 -5350];
        guesses.fwhm_guesses = [300 300 50];
        guesses.fwhm_lowerBounds = [0 0 0];
        guesses.fwhm_upperBounds = [1200 1200 1200];
        guesses.phase_guesses = [0 0 0];
        guesses.phase_lowerBounds = -inf*[1 1 1];
        guesses.phase_upperBounds = inf*[1 1 1];   
    end

    % scale areas for specific fid
    guesses.area_guess = (guesses.area_guess ./ max(guesses.area_guess)) * max(abs(fid));

    % time vector
    time = info.dwell*(0:info.nt-1);

    % intitialize fit
    fit = spec.NMR_TimeFit(fid,time,guesses.area_guess,...
        guesses.freq_guess,guesses.fwhm_guesses,guesses.phase_guesses,[],[]);
%     % set bounds
%     fit.setBounds(guesses.area_lowerBounds,guesses.area_upperBounds,...
%         guesses.freq_lowerBounds,guesses.freq_upperBounds,...
%         guesses.fwhm_lowerBounds,guesses.fwhm_upperBounds,...
%         guesses.phase_lowerBounds,guesses.phase_upperBounds);
    % fit
    fit.fitTimeDomainSignal();
    fit_result.fit = fit;

    % save initial and updated guesses
    fit_result.init_guesses = guesses;
    guesses.area_guess = fit.area; 
    guesses.freq_guess = fit.freq; 
    guesses.fwhm_guesses = fit.fwhm; 
    guesses.phase_guesses = fit.phase; 
    fit_result.new_guesses = guesses;

end

function images = decomp_fits(fits)
    %initialize
    % intensity
    gas = zeros(size(fits));
    membrane = zeros(size(fits));
    rbc = zeros(size(fits));
    % frequency
    gas_f = zeros(size(fits));
    membrane_f = zeros(size(fits));
    rbc_f = zeros(size(fits));
    % T2*
    gas_t2 = zeros(size(fits));
    membrane_t2 = zeros(size(fits));
    rbc_t2 = zeros(size(fits));
    
    % fill
    for x=1:size(gas,1)
        for y=1:size(gas,2)
            parfor z=1:size(gas,3)   
                % intensity
                gas(x,y,z) = fits{x,y,z}.fit.area(3) .* exp(1j*deg2rad(fits{x,y,z}.fit.phase(3)));
                membrane(x,y,z) = fits{x,y,z}.fit.area(2) .* exp(1j*deg2rad(fits{x,y,z}.fit.phase(2)));
                rbc(x,y,z) = fits{x,y,z}.fit.area(1) .* exp(1j*deg2rad(fits{x,y,z}.fit.phase(1)));
                % frequency
                gas_f(x,y,z) = fits{x,y,z}.fit.freq(3);
                membrane_f(x,y,z) = fits{x,y,z}.fit.freq(2);
                rbc_f(x,y,z) = fits{x,y,z}.fit.freq(1);
                % T2
                gas_t2(x,y,z) = 1/(pi * fits{x,y,z}.fit.fwhm(3))*1000; %in ms
                membrane_t2(x,y,z) = 1/(pi * fits{x,y,z}.fit.fwhm(2))*1000; %in ms
                rbc_t2(x,y,z) = 1/(pi * fits{x,y,z}.fit.fwhm(1))*1000; %in ms
                
            end
        end
    end

    % store
    % intensity
    images.gas = gas;
    images.membrane = membrane;
    images.rbc = rbc;
    % frequency
    images.gas_f = gas_f;
    images.membrane_f = membrane_f;
    images.rbc_f = rbc_f;
    % T2*
    images.gas_t2 = gas_t2;
    images.membrane_t2 = membrane_t2;
    images.rbc_t2 = rbc_t2;
end

function [maps, images] = calc_maps(images)
    % move quant images to maps since better fit
    fields = fieldnames(images);
    % assume proton image in first field and xenon images in second field
    sub_fields = fieldnames(images.(fields{2}));
    for im=1:length(sub_fields)
        if isreal(images.xe.(sub_fields{im})) %if real, then not complex image
            maps.(sub_fields{im}) = images.xe.(sub_fields{im});
            images.xe = rmfield(images.xe,sub_fields{im});
        end
    end

    % add in T2*/TE correction?
    fa_scale_factor = sind(20)/sind(0.5); % how to do accurately?
    maps.mem_quant = abs(images.xe.membrane ./ (images.xe.gas * fa_scale_factor));
    maps.rbc_quant = abs(images.xe.rbc ./ (images.xe.gas * fa_scale_factor));
    maps.rbc2mem_quant = abs(images.xe.rbc ./ images.xe.membrane);
end

function [images, maps] = interp_images(images, maps)
    %assume same fovs and H higher res
    [Xinterp, Yinterp, Zinterp] = meshgrid(...
        linspace(1,size(images.xe.gas,1),size(images.h,1)),... %x
        linspace(1,size(images.xe.gas,2),size(images.h,2)),... %y 
        linspace(1,size(images.xe.gas,3),size(images.h,3))); %z
    %interp xe images
    fields = fieldnames(images.xe);
    for im=1:length(fields)
        images.xe.(fields{im}) = interp3(images.xe.(fields{im}), Xinterp, Yinterp, Zinterp, 'makima');
    end
    %interp xe maps
    fields = fieldnames(maps);
    for im=1:length(fields)
        maps.(fields{im}) = interp3(maps.(fields{im}), Xinterp, Yinterp, Zinterp, 'makima');
    end


end

function h_image = h_recon(h_data, h_info)  
    pad_x = h_info.gen.kx_range(2)+h_info.gen.kx_range(1)+1;
    pad_y = h_info.gen.ky_range(2)+h_info.gen.ky_range(1)+1;
    pad_z = h_info.gen.kz_range(2)+h_info.gen.kz_range(1)+1;
    h_data = padarray(h_data,[pad_x pad_y pad_z],'pre');
    h_image = fftshift(fft(fftshift(h_data,1),[],1),1); % x
    h_image = fft(fftshift(h_image,2),[],2); % y
    h_image = fft(fftshift(h_image,3),[],3); % z
    h_image = squeeze(rssq(h_image, 4));
    out_size = round(size(h_image) ./ [h_info.gen.kx_oversample_factor h_info.gen.ky_oversample_factor h_info.gen.kz_oversample_factor]);
    crop_idx = centerCropWindow3d(size(h_image),out_size);
    h_image = imcrop3(h_image,crop_idx);
    % put in correct orientation for matlab with given acquisition
    h_image = flip(h_image,1);
end

function info = get_gen_list_info(list_file)
    % not set up for multiple mix/echo/loca
    if (strcmp(list_file(end-3:end),'data'))
        list_file(end-3:end) = 'list';
    end
    lines = readlines(list_file);
    lines(~startsWith(lines,'.    0    0    0  ')) = [];
    for str = 1:length(lines)
        % split at colon
        temp_str = split(lines(str),':');
        % make post-colon into array
        num = str2num(temp_str(2));
        % make pre-colon into var name
        temp_str = split(temp_str(1));
        temp_str = regexprep(temp_str,'-','_');
        temp_str = regexprep(temp_str,'0','z');
        temp_str = regexprep(temp_str,'1','o');
        temp_str = regexprep(temp_str,'2','t');
        % eval
        eval(strcat('info.',temp_str(end-1), ' = [', num2str(num),'];'))  
    end
end

function fig = make_fig(image_f,cmap_f,clims_f,alpha_f,image_b,cmap_b,clims_b)
    % set base values
    if size(image_f,4) ~= 3
        image_f = cat(4,image_f,image_f,image_f);
    end
    if size(image_b,4) ~= 3
        image_b = cat(4,image_b,image_b,image_b);
    end
    if isempty(cmap_f)
        cmap_f = 'hot';
    end
    if isempty(cmap_b)
        cmap_b = 'gray';
    end
    if isempty(clims_f)
        clims_f = [min(image_f(:)) max(image_f(:))];
    end
    if isempty(clims_b)
        clims_b = [min(image_b(:)) max(image_b(:))];
    end
    if isempty(alpha_f)
        alpha_f = 0.5*ones(size(image_f));
    end

    % normalize all to unit16
    bits = 16;
    bit_depth = 2^bits - 1;
    image_f = uint16(bit_depth * rescale(image_f,"InputMin",clims_f(1),"InputMax",clims_f(2)));
    image_b = uint16(bit_depth * rescale(image_b,"InputMin",clims_b(1),"InputMax",clims_b(2)));
    alpha_f(alpha_f>1) = 1;

    % convert to colormaps
    cmap_f = eval([cmap_f '(' num2str(bit_depth) ');']);
    cmap_b = eval([cmap_b '(' num2str(bit_depth) ');']);
    image_f_new = uint16(zeros(size(image_f)));
    image_b_new = uint16(zeros(size(image_b)));
    for x=1:size(image_f,1)
        for y=1:size(image_f,2)
            for z=1:size(image_f,3)
                image_f_new(x,y,z,:) = uint16(bit_depth*cmap_f(image_f(x,y,z,1)+1,:));
                image_b_new(x,y,z,:) = uint16(bit_depth*cmap_b(image_b(x,y,z,1)+1,:));
            end
        end
    end

    % linear comb
    fig = uint16(alpha_f.*double(image_f_new) + (1.-alpha_f).*double(image_b_new));

end

function save_images(images,path)
    fields = fieldnames(images);
    for f=1:length(fields)
        if isstruct(images.(fields{f}))
            sub_fields{f} = fieldnames(images.(fields{f}));
        end
    end
    for im=1:length(fields)
        if isempty(sub_fields{im})
            niftiwrite(images.(fields{im}), fullfile(path,['CSI_',fields{im}]))
        else
            parfor sub_im=1:length(sub_fields{im})
                if isreal(images.(fields{im}).(sub_fields{im}{sub_im}))
                    niftiwrite(images.(fields{im}).(sub_fields{im}{sub_im}), fullfile(path,['CSI_',fields{im},'_',sub_fields{im}{sub_im}]))
                else
                    niftiwrite(abs(images.(fields{im}).(sub_fields{im}{sub_im})), fullfile(path,['CSI_',fields{im},'_',sub_fields{im}{sub_im},'_mag']))
                    niftiwrite(angle(images.(fields{im}).(sub_fields{im}{sub_im})), fullfile(path,['CSI_',fields{im},'_',sub_fields{im}{sub_im},'_phase']))
                end
            end
        end
    end
end

function save_maps(maps,path)
    fields = fieldnames(maps);
    for im=1:length(fields)
        parfor sub_im=1:length(fields{im})
            niftiwrite(maps.(fields{im}), fullfile(path,['CSI_',fields{im},'_map']))
        end
    end
end

function save_figures(figures,path)
    options.color = true;
    options.compess = 'no';
    options.overwrite = true;
    options.message = false;

    fields = fieldnames(figures);
    for im=1:length(fields)
        saveastiff(permute(figures.(fields{im}),[1 2 4 3]),...
            fullfile(path,['CSI_', fields{im}, '_fig.tif']), options);
    end
end

function [results, maps] = quantify_results(maps, images)
    se = strel('cuboid', round(0.2*size(images.xe.gas)));
    % masks for calculations
    results.masks.gas = imerode(imbinarize(abs(images.xe.gas)),se);
    results.masks.mem = imerode(imbinarize(abs(images.xe.membrane)),se);
    results.masks.rbc = imerode(imbinarize(abs(images.xe.rbc)),se);
    results.masks.mem_uptake = (results.masks.gas & results.masks.mem);
    results.masks.rbc_transfer = (results.masks.gas & results.masks.rbc);
    results.masks.rbc2mem = (results.masks.rbc & results.masks.mem);

    % get stats
    results.stats.gas = datastats(abs(images.xe.gas(results.masks.gas)));
    results.stats.gas_f = datastats(maps.gas_f(results.masks.gas));
    results.stats.gas_t2 = datastats(maps.gas_t2(results.masks.gas));

    results.stats.mem_uptake = datastats(maps.mem_quant(results.masks.mem_uptake));
    results.stats.mem_f = datastats(maps.membrane_f(results.masks.mem));
    results.stats.mem_t2 = datastats(maps.membrane_t2(results.masks.mem));

    results.stats.rbc_transfer = datastats(maps.rbc_quant(results.masks.rbc_transfer));
    results.stats.rbc_f = datastats(maps.rbc_f(results.masks.rbc));
    results.stats.rbc_t2 = datastats(maps.rbc_t2(results.masks.rbc));

    results.stats.rbc2mem = datastats(maps.rbc2mem_quant(results.masks.rbc2mem));

    % alpha maps for figures
    max_alpha = 0.75;
    max_prct = 1;
    power = 1;
    width = 0.02;
    results.amap.gas = max_alpha*imgaussfilt3(rescale(abs(images.xe.gas),...
        "InputMin",0,"InputMax",prctile(abs(images.xe.gas(results.masks.gas)),max_prct)),...
        round(width*size(images.xe.gas))).^power;
    results.amap.mem = max_alpha*imgaussfilt3(rescale(abs(images.xe.membrane),...
        "InputMin",0,"InputMax",prctile(abs(images.xe.membrane(results.masks.mem)),max_prct)),...
        round(width*size(images.xe.gas))).^power;
    results.amap.rbc = max_alpha*imgaussfilt3(rescale(abs(images.xe.rbc),...
        "InputMin",0,"InputMax",prctile(abs(images.xe.rbc(results.masks.rbc)),max_prct)),...
        round(width*size(images.xe.gas))).^power;

end

% not implemented atm; only worry about if used more often
% function save_report(report,path)
% end
% 
% function rep = make_report(images, maps, figures, results)
%     rep = 1;
% end







