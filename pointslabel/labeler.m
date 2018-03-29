function [lab, anc] = labeler(lab, anc, w, h, nsweeps)
%LABELER Automatic text label placement, based on D3-Labeler
%
% [lab, anc] = labeler(lab, anc, w, h, nsweeps)
%
% This function optimizes the placement of text used to label points in a
% plot, minimizing the overlap between labels while keeping each label
% close to the point it is labeling, as well as constraining labels to stay
% within a specified bounding box.
%
% The code contained in this function is a translation of the the D3
% javascript plugin, D3-Labeler, written by Evan Wang
% (https://github.com/tinker10/D3-Labeler).  The input and output of this
% function is intended to mimic that version as closely as possible.  For a
% more Matlab-plot-friendly syntax, see the wrapper function pointslabel.m.
%
% Note that this function uses a matrix-style axis mode, where the
% coordinate system origin is in the upper left corner.
% 
% Input variables:
%
%   lab:        n x 1 labels structure, with the following fields
%               x:      initial x-coordinate of upper left corner of text
%                       label, in pixels  
%               y:      initial y-coordinate of upper left corner of text
%                       label, in pixels
%               width:  width of text label, in pixels
%               height: height of text label, in pixels
%               name:   text label
%
%   anc:        point anchors structure
%               x:      x-coordinate of center of anchor point, in pixels
%               y:      y-cooridnate of center of anchor point, in pixels
%               r:      radius of anchor point, in pixels
%
%   w:          width of canvas, in pixels
%
%   h:          height of canvas, in pixels
%
%   nsweeps:    number times to repeat random moves, per label
% 
% Output variables:
%
%   lab:        labels structure, with x- and y-coordinates updated with
%               the optimized positions.
%
%   anc:        point anchors structure

% Copyright 2017 Kelly Kearney


% Weights

w_len     = 0.2;  % leader line length 
w_inter   = 1.0;  % leader line intersection
w_lab2    = 30.0; % label-label overlap
w_lab_anc = 30.0; % label-anchor overlap
w_orient  = 3.0;  % orientation bias

max_move = 5.0;
max_angle = 0.5;
acc = 0;
rej = 0;

m = length(lab);
currT = 1.0;
initialT = 1.0;

for is = 1:nsweeps
    for im = 1:m 
        if rand(1) < 0.5
            mcmove(currT);
        else
            mcrotate(currT);
        end
    end
    currT = cooling_schedule(currT, initialT, nsweeps);
end

    function mcmove(currT) 
    % Monte Carlo translation move

        % select a random label
        i = randi([1 length(lab)], 1);

        % save old coordinates
        x_old = lab(i).x;
        y_old = lab(i).y;

        % old energy
        old_energy = energy(i);

        % random translation
        lab(i).x = lab(i).x + (rand(1) - 0.5) * max_move;
        lab(i).y = lab(i).y + (rand(1) - 0.5) * max_move;

        % hard wall boundaries
        if (lab(i).x > w) 
            lab(i).x = x_old;
        end
        if (lab(i).x < 0) 
            lab(i).x = x_old;
        end
        if (lab(i).y > h) 
            lab(i).y = y_old;
        end
        if (lab(i).y < 0) 
            lab(i).y = y_old;
        end

        % new energy
        new_energy = energy(i);

        % delta E
        delta_energy = new_energy - old_energy;

        if (rand(1) < exp(-delta_energy / currT))
            acc = acc + 1;
        else
            % move back to old coordinates
            lab(i).x = x_old;
            lab(i).y = y_old;
            rej = rej + 1;
        end
    end
    
    function mcrotate(currT) 
	% Monte Carlo rotation move

        % select a random label
        i = randi([1 length(lab)], 1);

        % save old coordinates
        x_old = lab(i).x;
        y_old = lab(i).y;

        % old energy
        old_energy = energy(i);

        % random angle
        angle = (rand(1) - 0.5) * max_angle;

        s = sin(angle);
        c = cos(angle);

        % translate label (relative to anchor at origin):
        lab(i).x = lab(i).x - anc(i).x;
        lab(i).y = lab(i).y - anc(i).y;

        % rotate label
        x_new = lab(i).x * c - lab(i).y * s;
        y_new = lab(i).x * s + lab(i).y * c;

        % translate label back
        lab(i).x = x_new + anc(i).x;
        lab(i).y = y_new + anc(i).y;

        % hard wall boundaries
        if (lab(i).x > w) 
            lab(i).x = x_old;
        end
        if (lab(i).x < 0) 
            lab(i).x = x_old;
        end
        if (lab(i).y > h) 
            lab(i).y = y_old;
        end
        if (lab(i).y < 0) 
            lab(i).y = y_old;
        end

        % new energy
        new_energy = energy(i);

        % delta E
        delta_energy = new_energy - old_energy;

        if (rand(1) < exp(-delta_energy / currT)) 
            acc = acc + 1;
        else
            % move back to old coordinates
            lab(i).x = x_old;
            lab(i).y = y_old;
            rej = rej + 1;
        end
    end

    function ener = energy(index)
    % energy function, tailored for label placement  

        ener = 0;
        dx = lab(index).x - anc(index).x;
        dy = anc(index).y - lab(index).y;
        dist = sqrt(dx * dx + dy * dy);
        overlap = true;
        amount = 0;
        theta = 0;

        % penalty for length of leader line

        if (dist > 0) 
            ener = ener + dist * w_len;
        end

        % label orientation bias

        dx = dx ./ dist;
        dy = dy ./ dist;
        if (dx > 0 && dy > 0) 
            ener = ener + 0 * w_orient; 
        elseif (dx < 0 && dy > 0) 
            ener = ener + 1 * w_orient; 
        elseif (dx < 0 && dy < 0) 
            ener = ener + 2 * w_orient;
        else
            ener = ener + 3 * w_orient;
        end

        x21 = lab(index).x;
        y21 = lab(index).y - lab(index).height + 2.0;
        x22 = lab(index).x + lab(index).width;
        y22 = lab(index).y + 2.0;

        % penalty for intersection of leader lines
        
        overlap = intersect(anc(index).x, lab(index).x, [anc.x], [lab.x], ...
                       anc(index).y, lab(index).y, [anc.y], [lab.y]);
        overlap(index) = 0; % skip self-overlap           
        ener = ener + w_inter .* sum(overlap);
        
        % penalty for label-label overlap
        
        x11 = [lab.x];
        y11 = [lab.y] - [lab.height] + 2.0;
        x12 = [lab.x] + [lab.width];
        y12 = [lab.y] + 2.0;
                   
        x_overlap = max(0, min(x12,x22) - max(x11,x21));
        y_overlap = max(0, min(y12,y22) - max(y11,y21));
        overlap_area = x_overlap .* y_overlap;
        overlap_area(index) = 0; % skip self-overlap
        
        ener = ener + sum(overlap_area) .* w_lab2;
        
%         for i = 1:m
%             if (i ~= index)
% 
%                 % penalty for intersection of leader lines
% 
% %                 overlap = intersect(anc(index).x, lab(index).x, anc(i).x, lab(i).x,...
% %                                     anc(index).y, lab(index).y, anc(i).y, lab(i).y);
%                 if (overlap(i)) 
%                     ener = ener + w_inter;
%                 end
% 
%                 % penalty for label-label overlap
% 
%                 x11 = lab(i).x;
%                 y11 = lab(i).y - lab(i).height + 2.0;
%                 x12 = lab(i).x + lab(i).width;
%                 y12 = lab(i).y + 2.0;
% 
%                 x_overlap = max(0, min(x12,x22) - max(x11,x21));
%                 y_overlap = max(0, min(y12,y22) - max(y11,y21));
%                 overlap_area = x_overlap * y_overlap;
%                 ener = ener + (overlap_area * w_lab2);
%             end
% 
%             % penalty for label-anchor overlap
% 
%             x11 = anc(i).x - anc(i).r;
%             y11 = anc(i).y - anc(i).r;
%             x12 = anc(i).x + anc(i).r;
%             y12 = anc(i).y + anc(i).r;
%             x_overlap = max(0, min(x12,x22) - max(x11,x21));
%             y_overlap = max(0, min(y12,y22) - max(y11,y21));
%             overlap_area = x_overlap * y_overlap;
%             ener = ener + (overlap_area * w_lab_anc);
%         end
    end

end


function inter = intersect(x1, x2, x3, x4, y1, y2, y3, y4) 
% returns true if two lines intersect, else false
% from http:%paulbourke.net/geometry/lineline2d/

    denom = (y4 - y3) .* (x2 - x1) - (x4 - x3) .* (y2 - y1);
    numera = (x4 - x3) .* (y1 - y3) - (y4 - y3) .* (x1 - x3);
    numerb = (x2 - x1) .* (y1 - y3) - (y2 - y1) .* (x1 - x3);

    % Is the intersection along the the segments 
    mua = numera ./ denom;
    mub = numerb ./ denom;

    inter = ~(mua < 0 | mua > 1 | mub < 0 | mub > 1);

end

function t = cooling_schedule(currT, initialT, nsweeps) 
% Linear cooling schedule
    t = (currT - (initialT / nsweeps));
end



