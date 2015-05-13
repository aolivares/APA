function [xmin, fmin, ct] = ANMS(myfunction, xinit, tol, max_feval)
%
% FUNCTION ANMS uses the Adaptive Nelder-Mead Simplex Algorithm to minimize
% 'myfunction'. This algorithm uses the adaptive parameters introduced in
% Fuchang Gao and Lixing Han "Implementing the Nelder-Mead simplex 
% algorithm with adaptive parameters" Computational Optimization and 
% Applications, appeared online May 4, 2010. It also appliesa a large 
% initial simplex. 
%  
% - INPUT PARAMETERS:
%   |_ 'myfunction': Objective function to be minimized.
%   |_ 'xinit': Initial guess.
%   |_ 'tol': Tolerance for termination (Recommended value: 10^-4);
%   |_ 'max_feval': Maximum number of function evaluations. (Remark: In 
%                   myfunction(x), x must be a row vector).
%
% - OUTPUT PARAMETERS: 
%   |_ 'xmin': Approximate optimal solution at termination.
%   |_ 'fmin': Minimum function value at termination.
%   |_ 'ct': Number of function evaluations at termination.

% -------------------------------------------------------------------------
% |||||||||||||||||||||| COPYRIGHT AND AUTHORSHIP |||||||||||||||||||||||||
% -------------------------------------------------------------------------
% *************************************************************************
% AUTHORS:
%   |_ Fuchang Gao* and Lixing Han** (August, 2010), Original version.
%       * Entity: Department of Mathematics, University of Idaho.
%       * Contact: fuchang@uidaho.edu
%       ** Entity: Department of Mathematics, University of Michigan.
%       ** Contact: lxhan@umflint.edu
%   |_ Dr. Alberto Olivares (Reorganization of code):
%       * Entity:   Department of Signal Theory, Telematics and
%                   Communications, University of Granada, Spain.
%       * Contact:  aolivares@ugr.es
%
% LAST MODIFICATION: 
%   |_ Date: 01/22/2015 (mm/dd/yy). 
%   |_ Location: Research Centre for Information and Communications 
%                Technologies of the University of Granada, Spain. 
%                (CITIC-UGR).
% *************************************************************************
%
%  Copyright (c) 2015, Fuchang Gao and Lixing Han and Alberto Olivares. 
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the 
%       distribution.
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
%  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
%  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
%  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
%  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
%  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
%  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
%  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------

% x0 is a row vector.
x0 = xinit(:)';  
myfunction = fcnchk(myfunction);

% Dimension of the problem.
dim = max(size(x0)); 

% Set up adaptive parameters.
alpha = 1; 
beta = 1 + 2 / dim; 
gamma = 0.75 - 0.5 / dim; 
delta = 1 - 1 / dim;

% Construct the initial simplex: Large initial simplex is used
scalefactor = min(max(max(abs(x0)), 1), 10);
D0 = eye(dim);
D0(dim + 1, :) = (1 - sqrt(dim + 1)) / dim * ones(1, dim);
for i = 1 : dim + 1
    X(i, :) = x0 + scalefactor * D0(i, :);
    FX(i) = feval(myfunction, X(i, :));
end;
ct = dim + 1;
[FX, I] = sort(FX);
X = X(I, :);

% Main iteration.
while max(max(abs(X(2 : dim + 1, :) - X(1 : dim, :)))) >= scalefactor * tol 
    if ct > max_feval
        break;
    end   
    % Centroid of the dim best vertices.
    M = mean(X(1 :dim, :));  
    xref = (1 + alpha) * M - alpha * X(dim + 1, :);
    Fref = feval(myfunction, xref);
    ct = ct + 1;
    if Fref < FX(1)
        % expansion
        xexp = (1 + alpha * beta) * M - alpha * beta * X(dim + 1, :);
        Fexp = feval(myfunction, xexp);
        ct = ct + 1;
        if Fexp < Fref
            X(dim + 1, :) = xexp;
            FX(dim + 1) = Fexp;
        else
            X(dim + 1, :)= xref;
            FX(dim + 1) = Fref;
        end;
    else
        if Fref < FX(dim)
            % accept reflection point.
            X(dim + 1, :) = xref;
            FX(dim + 1) = Fref;
        else 
            if Fref < FX(dim + 1)
                % Outside contraction
                xoc = (1 + alpha * gamma) * M - alpha * gamma * ...
                    X(dim + 1, :);
                Foc = feval(myfunction, xoc);
                ct = ct + 1;           
                if Foc <= Fref
                    X(dim + 1, :) = xoc;
                    FX(dim + 1) = Foc;
                else
                    % shrink
                    for i = 2 : dim + 1
                        X(i, :) = X(1, :) + delta * (X(i, :) - X(1, :));
                        FX(i) = feval(myfunction, X(i, :));
                    end;
                    ct = ct + dim;
                end;
            else
                %inside contraction
                xic = (1 - gamma) * M + gamma * X(dim + 1, :);
                Fic = feval(myfunction, xic);
                ct = ct + 1;
                if Fic < FX(dim + 1)
                    X(dim + 1, :) = xic;
                    FX(dim + 1) = Fic;
                else
                    % shrink
                    for i = 2 : dim + 1
                        X(i, :) = X(1, :) + delta * (X(i, :) - X(1, :));
                        FX(i) = feval(myfunction, X(i, :));
                    end;
                    ct = ct + dim;
                end;
            end;
        end;
    end;
    [FX, I] = sort(FX);
    X = X(I, :);
end
xmin = X(1 ,:);
fmin = FX(1);
