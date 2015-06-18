% This file is part of the SPLINTER library.
% Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
%
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

classdef RadialBasisFunctionType
    enumeration
       Thin_plate_spline,
       Multiquadric,
       Inverse_quadric,
       Inverse_multiquadric,
       Gaussian
    end
end