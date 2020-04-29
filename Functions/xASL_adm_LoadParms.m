function [Parms,x,Oldx] = xASL_adm_LoadParms(ParmsPath, x, bVerbose)
%xASL_adm_LoadParmsMat Loads parameters from .mat parameter file
% NB: future version should include JSON
% If x are provided, this function will compare x with ParmsPath

% strMatError contains errors about non-existent mat. Only report this error when JSON is also missing


%% ------------------------------------------------------------------------
%% 0) Admin
Parms = struct; % default

if nargin<1 || isempty(ParmsPath)
    error('ParmsPath was not specified');
end
if nargin<2 || isempty(x)
    warning('No x struct fields available, can be missing some information for quantification');
    x = struct; % dummy variable
end
if nargin<3 || isempty(bVerbose)
    bVerbose = true;
end

strMatError = ''; % initializing string containing errors/warning to combine & throw as warning at end

[Fpath, Ffile, Fext] = fileparts(ParmsPath);

%% ------------------------------------------------------------------------
%% 1) Load .mat parameter file (if exists)
if ~exist(ParmsPath, 'file') && isfield(x.Q,'ASL')
    warning('No ParmsPath found, trying to use general parameters!');
    
    if ~isempty(regexp(Ffile,'ASL4D'))
        if bVerbose; fprintf('Use generic ASL parameters\n'); end
        ImType = 'ASL';
    elseif ~isempty(regexp(Ffile,'M0'))
        if bVerbose; fprintf('Use generic M0 parameters\n'); end
        ImType = 'M0';
    end
    
    FieldsGeneric = fields(x.Q.(ImType));
    for iA=1:length(FieldsGeneric)
        if ~isfield(x.Q,FieldsGeneric{iA})
            x.Q.(FieldsGeneric{iA}) = x.Q.(ImType).(FieldsGeneric{iA});
        end
    end
    
elseif ~exist(ParmsPath, 'file')
	strMatError = '_parms.mat file missing';
    
elseif strcmp(Fext,'.mat')
    Parms = load(ParmsPath,'-mat');
    if  isfield(Parms,'parms')
        Parms = Parms.parms;
    else
		strMatError = sprintf('ERROR in module_ASL: could not load M0 parameters from:\n%s\n',M0_parms_file);
    end
end

%% ------------------------------------------------------------------------
%% 2) Load JSON file (if exists)
% Define JSON path
if ~isempty(regexp(Ffile,'ASL4D'))
    JSONPath = fullfile(Fpath, 'ASL4D.json');
elseif ~isempty(regexp(Ffile,'M0'))
    JSONPath = fullfile(Fpath, 'M0.json');
elseif ~isempty(regexp(Ffile,'func'))
    JSONPath = xASL_adm_GetFileList(Fpath,'func.*bold\.json',  'FPList', [0 Inf]);
    if ~isempty(JSONPath)
        JSONPath = JSONPath{1};
    end
elseif ~isempty(regexp(Ffile,'dwi'))
    JSONPath = xASL_adm_GetFileList(Fpath,'dwi.*dwi\.json',  'FPList', [0 Inf]);
    if ~isempty(JSONPath)
        JSONPath = JSONPath{1};
    end    
else
    warning('Could not define JSON path');
end

% Load JSON file
VendorRestore = false; % default
if exist(JSONPath,'file') % According to the BIDS inheritance principle, the JSON values overwrite the existing values
    % But this can have some exceptions, e.g. the specification of sequence
    % in x.Vendor
    
    if isfield(x,'Vendor') && (strcmp(x.Vendor,'GE_product') || strcmp(x.Vendor,'GE_WIP'))
        VendorBackup = x.Vendor;
        VendorRestore = true;
    end
    
	JSONParms = spm_jsonread(JSONPath);
    
    ParmsFields = {'EchoTime' 'RepetitionTime' 'RescaleSlopeOriginal' 'RescaleSlope'    'MRScaleSlope'...
        'AcquisitionTime' 'RescaleIntercept'        'readout_dim'       'Vendor'        'readout_dim'...
        'Vendor' 'BackGrSupprPulses' 'LabelingType' 'Initial_PLD' 'LabelingDuration' 'SliceReadoutTime' 'M0'};
    JSONFields  = {'EchoTime' 'RepetitionTime' 'PhilipsRescaleSlope'  'PhilipsRWVSlope' 'PhilipsScaleSlope'...
        'AcquisitionTime' 'PhilipsRescaleIntercept' 'MRAcquisitionType' 'Manufacturer'  'readout_dim'...
        'Vendor' 'BackGrSupprPulses' 'LabelingType' 'Initial_PLD' 'LabelingDuration' 'SliceReadoutTime' 'M0'};
    
    for iP=1:length(JSONFields)
        if isfield(JSONParms,JSONFields{iP})
            Parms.(ParmsFields{iP}) = JSONParms.(JSONFields{iP});
            if ~isempty(regexp(JSONFields{iP}, 'EchoTime'))  % JSON is in seconds, DICOM in ms
                if Parms.(ParmsFields{iP})>0.001 && Parms.(ParmsFields{iP})<1
                    % we expect 1 < TE(ms) < 1000
                    Parms.(ParmsFields{iP}) = Parms.(ParmsFields{iP})*1000;
                else
                    warning('JSON EchoTime was not according to BIDS (seconds), check quantification');
                    fprintf('DICOM uses TE & TR in ms, BIDS JSON in s');
                end
            elseif ~isempty(regexp(JSONFields{iP}, 'RepetitionTime'))  % JSON is in seconds, DICOM in ms
                if Parms.(ParmsFields{iP})>0.1 && Parms.(ParmsFields{iP})<100
                    % we expect 100 < TR(ms) < 100,000
                    Parms.(ParmsFields{iP}) = Parms.(ParmsFields{iP})*1000;
                else
                    warning('JSON EchoTime was not according to BIDS (seconds), check quantification');
                    fprintf('DICOM uses TE & TR in ms, BIDS JSON in s');
                end                    
            elseif strcmp(ParmsFields{iP},'AcquisitionTime')
                tP = xASL_adm_CorrectName(Parms.AcquisitionTime, 2);
                tP = xASL_adm_ConvertNr2Time(xASL_adm_ConvertTime2Nr(tP));
                Parms.AcquisitionTime = [tP(1:6) '.' tP(7:8)];
            end
        end
    end
else
    strMatError = [strMatError ', JSON file missing'];
end
    
%% ------------------------------------------------------------------------
%% 3) Deal with warnings
if ~exist('Parms','var') && isempty(strMatError)
    Parms = struct; 
    warning('parms seem missing, something wrong with parmsfile?');
elseif ~exist('Parms','var') && ~isempty(strMatError)
    warning([strMatError ', skipping...']);
    return;
elseif ~isempty(strMatError)
    if bVerbose; fprintf('%s\n',[strMatError ', but we try using the other source (*_parms.mat/*.json)']); end
end

if VendorRestore
    x.Vendor = VendorBackup;
    Parms.Vendor = VendorBackup;
end

%% ------------------------------------------------------------------------
%% 4) Fix M0 parameter if not set
if ~isfield(x,'M0')
    if xASL_exist(fullfile(Fpath, 'M0.nii'),'file') && (exist(fullfile(Fpath, 'M0.json'),'file') || exist(fullfile(Fpath, 'M0_parms.mat'),'file') )
        x.M0 = 'separate_scan';
        if bVerbose; fprintf('%s\n',['M0 parameter was missing, set to ' x.M0]); end
    elseif isfield(Parms,'BackGrSupprPulses') && Parms.BackGrSupprPulses==0
        x.M0 = 'UseControlAsM0';
        if bVerbose; fprintf('%s\n',['M0 parameter was missing, set to ' x.M0]); end
    else
        if bVerbose; fprintf('%s\n','M0 parameter was missing, OR didnt find M0 scan, AND BackGrSupprPulses wasnt set to 0...'); end
    
    end
end

if ~exist('Parms','var')
    Parms = struct; 
    warning('parms seem missing, something wrong with parmsfile?');
end

%% ------------------------------------------------------------------------
%% 5) Check erroneous scale slope 
if  isfield(Parms,'RescaleSlope') && isfield(Parms,'RescaleSlopeOriginal')
    RelDiff = abs(100* (Parms.RescaleSlopeOriginal-Parms.RescaleSlope) / ((Parms.RescaleSlopeOriginal+Parms.RescaleSlope)*0.5));
    if  RelDiff>0.05
        warning(['RescaleSlope & RescaleSlopeOriginal were not identical in ' ParmsPath]);
        fprintf('%s\n', ', check if scaling was applied (especially on Philips scanner enhanced DICOMs');
        fprintf('%s\n',['Difference was ' xASL_num2str(RelDiff) '%']);
    end
end


%% ------------------------------------------------------------------------
%% 6) Backwards compatibility section
% Input all fields from the Parms into the x structure, backup those that were already existing (inheritance principle)
Oldx            = struct;
ParmsFieldNames = fieldnames(Parms);
xFieldNames     = fieldnames(x);

for iP=1:length(ParmsFieldNames)
    % Check first if already exists, to backup
    for iS=1:length(xFieldNames)
        if  strcmp(ParmsFieldNames{iP},xFieldNames{iS})
            Oldx.(xFieldNames{iS}) = x.(xFieldNames{iS}); % backup
        end
    end
    x.(ParmsFieldNames{iP}) = Parms.(ParmsFieldNames{iP});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % QUICK & DIRTY FIX
    x.Q.(ParmsFieldNames{iP}) = Parms.(ParmsFieldNames{iP});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Move quantification parameters to the Q (quantification) subfield, for
% backward compatibility
Qfields = {'BackGrSupprPulses' 'LabelingType' 'Initial_PLD' 'LabelingDuration' 'SliceReadoutTime' 'Lambda' 'T2art' 'BloodT1' 'TissueT1' 'nCompartments'};
for iField=1:length(Qfields)
    if isfield(x,Qfields{iField})
        if isfield(x.Q,(Qfields{iField})) && ~strcmp(x.Q.(Qfields{iField}), x.(Qfields{iField}))
            warning(['Overwriting x.Q.' Qfields{iField} '=' x.Q.(Qfields{iField}) ', with x.' Qfields{iField} '=' x.(Qfields{iField})]);
        end
        
        x.Q.(Qfields{iField}) = x.(Qfields{iField});
        x = rmfield(x, Qfields{iField});
    end
end



end
