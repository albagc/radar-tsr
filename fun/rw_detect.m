function [Mrw, mrw_spe, mrw_t2, Mcw] = rw_detect(spe, t2, ucl_spe, ucl_t2, fac_spe, fac_t2, alpha, e, M, mtrue, Mtrue, dplot)
%--------------------------------------------------------------------------
% rw_detect - Detection of rowwise outliers based on SPE and T2 statistics
%--------------------------------------------------------------------------
%
% Inputs:
%   spe: SPE statistics (vector)
%   t2: T2 statistics (vector)
%   ucl_spe: Upper control limit for SPE (scalar)
%   ucl_t2: Upper control limit for T2 (scalar)
%   fac_spe: Factor for SPE cutoff calculation (default: 2.5)
%   fac_t2: Factor for T2 cutoff calculation (default: 2.5)
%   alpha: False positive rate (default: 0.01)
%   e: Matrix indicating cellwise outliers (optional)
%   M: Matrix indicating missing values (optional)
%   mtrue: Matrix indicating true rowwise outliers (optional)
%   Mtrue: Matrix indicating true missing values (optional)
%   dplot: Logical flag for diagnostic plots (default: false)
%
% Outputs:
%   Mrw: Logical vector indicating rowwise outliers
%   mrw_spe: Logical vector indicating rowwise outliers based on SPE
%   mrw_t2: Logical vector indicating rowwise outliers based on T2
%   Mcw: Logical matrix indicating cellwise outliers

arguments
    spe
    t2
    ucl_spe
    ucl_t2
    fac_spe = 2.5
    fac_t2 = 2.5
    alpha = 0.01
    e = nan
    M = false(size(e))
    mtrue = false(size(spe))
    Mtrue = false(size(M))
    dplot = false
end

% Number of observations expected above the UCL
n_falserate = round(alpha * numel(spe));

% Check if the number of observations above the UCL exceeds the false positive rate
if sum(spe > ucl_spe) > n_falserate
    Mrw_spe_cutoff = spe > (fac_spe * ucl_spe);
else
    Mrw_spe_cutoff = false(size(spe));
end

if sum(t2 > ucl_t2) > n_falserate
    Mrw_t2_cutoff = t2 > (fac_t2 * ucl_t2);
else
    Mrw_t2_cutoff = false(size(t2));
end

Mrw_cutoff = Mrw_t2_cutoff | Mrw_spe_cutoff;

if ~isempty(e) && ~any(islogical(e))
    Mcw = univarfilter(e, alpha, 'rob', zeros(size(e)), 'auto', false, true);
    Mcw(M) = false;
elseif islogical(e)
    Mcw = e;
end

Mrw_alpha_ucl_spe = spe > ucl_spe;
ind_alpha_ucl_spe = find(Mrw_alpha_ucl_spe);
[~, i_sort_spe] = sort(spe(Mrw_alpha_ucl_spe), 'ascend');

if numel(i_sort_spe) <= n_falserate
    mrw_spe = Mrw_spe_cutoff;
else
    Mrw_alpha_ucl_spe(ind_alpha_ucl_spe(i_sort_spe(1:n_falserate))) = false;
    
    if fac_spe == 1
        mrw_spe = Mrw_spe_cutoff | Mrw_alpha_ucl_spe;
    elseif fac_spe == 2.5
        mrw_spe = Mrw_spe_cutoff;
    end
end

Mrw_alpha_ucl_t2 = t2 > ucl_t2;
ind_alpha_ucl_t2 = find(Mrw_alpha_ucl_t2);
[~, i_sort_t2] = sort(t2(Mrw_alpha_ucl_t2), 'ascend');

if numel(i_sort_t2) <= n_falserate
    mrw_t2 = Mrw_t2_cutoff;
else
    Mrw_alpha_ucl_t2(ind_alpha_ucl_t2(i_sort_t2(1:n_falserate))) = false;
    
    if fac_t2 == 1
        mrw_t2 = Mrw_t2_cutoff | Mrw_alpha_ucl_t2;
    elseif fac_t2 == 2.5
        mrw_t2 = Mrw_t2_cutoff;
    end
end

if numel(e) > 1
    [m_ecw, rw_efilt] = rowfiltcw(Mcw, false(size(Mcw)), alpha, 'rob');
    m_ecw = mrw_spe & m_ecw;
    m_cw = mrw_spe & ~Mrw_cutoff;
    Mcw(~m_cw, :) = false;
    Mrw = (mrw_spe | mrw_t2) & ~m_cw;
else
    Mcw = nan;
    Mrw = mrw_spe | mrw_t2;
end

end
