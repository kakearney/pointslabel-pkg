function [ht, ha] = pointslabel(xanchor, yanchor, labels, varargin)
%POINTSLABEL Add text labels to points, using optimized positioning
%
% [ht, ha] = pointslabel(xanchor, yanchor, labels, varargin)
%
% This function adds smart text labels to points, attempting to minimize
% overlap between labels while still keeping each text label close to the
% point it is labeling.
%
% Input variables:
%
%   xanchor:    n x 1 array, x-coordinates of the text anchors, i.e. the
%               points to be labeled
%
%   yanchor:    n x 1 array, y-coordinates of the text anchors
%
%   labels:     n x 1 array of text labels.  Each label can be either a
%               string, a character array, or a cell array of character
%               arrays (for multi-line labels).
%
% Optional input variables (passed as parameter/value pairs):
%
%   nsweeps:    number of sweeps used by the text placement optimization
%               alorithm [1000]
%
%   bbox:       area to use as a bounding box for the text ['axis']
%               'axis':     Text labels will be constrained to stay within
%                           the specified axis 
%
%               'figure':   Text labels will be constrained to stay within
%                           the parent figure of the specified axis (WIP)
%
%               [l r w h]:  Text labels will be contrained within the
%                           specified bounding box, in normalized figure
%                           units (TODO: not yet coded)
%
%   axis:       handle of axis where text labels will be added [gca]
%
%   radius:     radius (in pixels) of the anchor points.  Text labels will
%               remain at minimum this distance from each point.

% Copyright 2017 Kelly Kearney

istxt = @(x,n) ischar(x) || isstring(x) || (iscellstr(x) && numel(x)==n);

% Check x and y labels

validateattributes(xanchor, {'numeric'}, {});
npt = numel(xanchor);
validateattributes(yanchor, {'numeric'}, {'numel', npt});

xanchor = xanchor(:);
yanchor = yanchor(:);

% Parse optional

p = inputParser;
p.addParameter('nsweeps',     1000,      @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
p.addParameter('bbox',        'axis',    @(x) validateattributes(x, {'char', 'numeric'}));
p.addParameter('axis',        gca);
p.addParameter('radius', 0);


textProps = {'FontName','FontSize','FontUnits','FontWeight','Interpreter'};
tprop = cellfun(@(x) get(0, ['DefaultText' x]), textProps, 'uni', 0);

for ii = 1:length(textProps)
    p.addParameter(textProps{ii}, tprop{ii});
end
p.parse(varargin{:});
Opt = p.Results;

for ii = 1:length(textProps)
    if iscell(Opt.(textProps{ii}))
        if ~(numel(Opt.(textProps{ii})) == npt)
            error('Text properties must be scalar or cell array with same size as x and y');
        end
    else
        Opt.(textProps{ii}) = repmat({Opt.(textProps{ii})}, npt, 1);
    end
end

% TODO: Check that all labels are strings/character arrays or cell arrays of strings

if npt == 1 && ~(iscell(labels) && numel(labels) == 1)
    labels = {labels};
end
labels = labels(:);

% Calculate text position and extent and canvas width/height in pixels

hfig = ancestor(Opt.axis, 'figure');
if ischar(Opt.bbox)
    switch Opt.bbox
        
        %--------------------
        % Canvas is plot axis
        %--------------------
        
        case 'axis'
            
            % Convert to pixel coordinates, using image-style coordinates
            % with origin at axis upper left
            
            units = get(Opt.axis, 'units');
            set(Opt.axis, 'units', 'pixels');
            axpos = plotboxpos(Opt.axis);
            set(Opt.axis, 'units', units);
            
            xlim = get(Opt.axis, 'xlim');
            ylim = get(Opt.axis, 'ylim');
            
            w = axpos(3);
            h = axpos(4);
            
            xfwd = @(xtxt) (xtxt - xlim(1))./(diff(xlim)) * w; % Axis coordinates to ij-pixel coordinate
            yfwd = @(ytxt) (ylim(2) - ytxt)./(diff(ylim)) * h;
            xrev = @(xpix) (diff(xlim)*xpix)./w + xlim(1); % ij-pixel coords to axis coords
            yrev = @(ypix) ylim(2) - (diff(ylim)*ypix)./h;
            
            x = xfwd(xanchor);
            y = yfwd(yanchor);
            
%             x = ((xanchor - xlim(1))./diff(xlim)) * w;
%             y = ((ylim(2) - yanchor)./diff(ylim)) * h;
            
            % Calculate label extent, in pixels
            
            [lw, lh] = deal(nan(npt,1));
            
            pv = cellfun(@(x) Opt.(x)(:), textProps, 'uni', 0);
            pv = cat(2, pv{:});
            
            for ii = 1:npt
                htxt = text(xlim(1), ylim(1), labels{ii}, 'visible', 'off');
                set(htxt, textProps', pv(ii,:));
                set(htxt, 'units', 'pixels');
                ex = get(htxt, 'extent');
                lw(ii) = ex(3);
                lh(ii) = ex(4);
                delete(htxt);
            end
            
            % Call labeler
            
            anc = struct('x', num2cell(x), 'y', num2cell(y), 'r', Opt.radius); 
            lab = struct('x', num2cell(x), 'y', num2cell(y), ...
                'width', num2cell(lw), 'height', num2cell(lh), 'name', labels);
                      
            [lab, anc] = labeler(lab, anc, w, h, Opt.nsweeps);
            
            % Calculate leader line intersection point
            
            for ii = 1:npt
                
                x1 = lab(ii).x + [0 0 1 1 0].*lab(ii).width;
                y1 = lab(ii).y + [0 1 1 0 0].*lab(ii).height;
                
                xy = distance2curve([x1;y1]', [anc(ii).x anc(ii).y]);
                anc(ii).xlead = xy(1);
                anc(ii).ylead = xy(2);
                
%                 x2 = [lab(ii).x+lab(ii).width/2 anc(ii).x];
%                 y2 = [lab(ii).y+lab(ii).height/2 anc(ii).y];
%                 [x3,y3] = intersections([x1 NaN x2], [y1 NaN y2]);
%                 
%                 xy = setdiff([x3 y3], [x1' y1'], 'rows');
%                 if isempty(xy)
%                     anc(ii).xlead = NaN;
%                     anc(ii).ylead = NaN;
%                 else
%                     anc(ii).xlead = xy(1);
%                     anc(ii).ylead = xy(2);
%                 end
                    
            end
            
            % Translate back to data units
            
            xtxt  = xrev([lab.x]);
            xanc  = xrev([anc.x]);
            xlead = xrev([anc.xlead]);
            
            ytxt  = yrev([lab.y]);
            yanc  = yrev([anc.y]);
            ylead = yrev([anc.ylead]);
            
            
            
%             xtxt = ([lab.x]./w).*diff(xlim) + xlim(1);
%             xanc = ([anc.x]./w).*diff(xlim) + xlim(1);
%             xlead = ([anc.xlead]./w).*diff(xlim) + xlim(1);
%    
%             ytxt = ([lab.y]./h).*diff(ylim) + ylim(1);
%             yanc = ([anc.y]./h).*diff(ylim) + ylim(1);
%             ylead = ([anc.ylead]./h).*diff(ylim) + ylim(1);
            

            ht = text(xtxt, ytxt, labels, 'horiz', 'left', 'vert', 'top');
            set(ht, textProps', pv); 
            ha = line([xanc; xlead], [yanc; ylead], 'color', [0.8 0.8 0.8]);
            
        case 'figure' % Canvas is figure
            
            hfig = ancestor(Opt.axis, 'figure');
           
            figpos = getpos(hfig, 'px');
            w = figpos(3);
            h = figpos(4);
            
            figunit = get(hfig, 'Units');
            set(hfig, 'units', 'pixels');
            [x, y] = axescoord2figurecoord(xanchor, yanchor, Opt.axis);
            
            error('Haven''t finished figure reverse yet');
           
            
        otherwise % Custom canvas (in plot axis units)
            error('bbox must be ''axis'', ''figure'', or a 1x4 bounding box array');
            
    end
else
    % Coordinates are relative to region in axis coordinates
    validateattributes(Opt.box, {'numeric'}, {'positive', 'size', [1 4]});
    
    error('Haven''t coded custom canvas yet');
end





% %%%%%%%%%%%%%%%%%%
% 
% function lab, anc = labeler(lab, anc, w, h, nsweeps)
% % Main simulated annealing
% 
%     % Weights
% 
%     w_len     = 0.2;  % leader line length 
%     w_inter   = 1.0;  % leader line intersection
%     w_lab2    = 30.0; % label-label overlap
%     w_lab_anc = 30.0; % label-anchor overlap
%     w_orient  = 3.0;  % orientation bias
% 
%     max_move = 5.0;
%     max_angle = 0.5;
%     acc = 0;
%     rej = 0;
% 
%     m = length(lab);
%     currT = 1.0;
%     initialT = 1.0;
% 
%     for is = 1:nsweeps
%         for im = 1:m 
%             if rand(1) < 0.5
%                 mcmove(currT);
%             else
%                 mcrotate(currT);
%             end
%         end
%         currT = cooling_schedule(currT, initialT, nsweeps);
%     end
% 
%     function mcmove(currT) 
%     % Monte Carlo translation move
% 
%         % select a random label
%         i = randi([1 length(lab)], 1);
% 
%         % save old coordinates
%         x_old = lab(i).x;
%         y_old = lab(i).y;
% 
%         % old energy
%         old_energy = energy(i);
% 
%         % random translation
%         lab(i).x = lab(i).x + (rand(1) - 0.5) * max_move;
%         lab(i).y = lab(i).y + (rand(1) - 0.5) * max_move;
% 
%         % hard wall boundaries
%         if (lab(i).x > w) 
%             lab(i).x = x_old;
%         end
%         if (lab(i).x < 0) 
%             lab(i).x = x_old;
%         end
%         if (lab(i).y > h) 
%             lab(i).y = y_old;
%         end
%         if (lab(i).y < 0) 
%             lab(i).y = y_old;
%         end
% 
%         % new energy
%         new_energy = energy(i);
% 
%         % delta E
%         delta_energy = new_energy - old_energy;
% 
%         if (rand(1) < exp(-delta_energy / currT))
%             acc = acc + 1;
%         else
%             % move back to old coordinates
%             lab(i).x = x_old;
%             lab(i).y = y_old;
%             rej = rej + 1;
%         end
%     end
%     
%     function mcrotate(currT) 
% 	% Monte Carlo rotation move
% 
%         % select a random label
%         i = randi([1 length(lab)], 1);
% 
%         % save old coordinates
%         x_old = lab(i).x;
%         y_old = lab(i).y;
% 
%         % old energy
%         old_energy = energy(i);
% 
%         % random angle
%         angle = (rand(1) - 0.5) * max_angle;
% 
%         s = sin(angle);
%         c = cos(angle);
% 
%         % translate label (relative to anchor at origin):
%         lab(i).x = lab(i).x - anc(i).x;
%         lab(i).y = lab(i).y - anc(i).y;
% 
%         % rotate label
%         x_new = lab(i).x * c - lab(i).y * s,
%         y_new = lab(i).x * s + lab(i).y * c;
% 
%         % translate label back
%         lab(i).x = x_new + anc(i).x;
%         lab(i).y = y_new + anc(i).y;
% 
%         % hard wall boundaries
%         if (lab(i).x > w) 
%             lab(i).x = x_old;
%         end
%         if (lab(i).x < 0) 
%             lab(i).x = x_old;
%         end
%         if (lab(i).y > h) 
%             lab(i).y = y_old;
%         end
%         if (lab(i).y < 0) 
%             lab(i).y = y_old;
%         end
% 
%         % new energy
%         new_energy = energy(i);
% 
%         % delta E
%         delta_energy = new_energy - old_energy;
% 
%         if (rand(1) < Math.exp(-delta_energy / currT)) 
%             acc = acc + 1;
%         else
%             % move back to old coordinates
%             lab(i).x = x_old;
%             lab(i).y = y_old;
%             rej = rej + 1;
%         end
%     end
% 
%     function ener = energy(index)
%     % energy function, tailored for label placement  
% 
%         ener = 0;
%         dx = lab(index).x - anc(index).x;
%         dy = anc(index).y - lab(index).y;
%         dist = sqrt(dx * dx + dy * dy);
%         overlap = true;
%         amount = 0;
%         theta = 0;
% 
%         % penalty for length of leader line
% 
%         if (dist > 0) 
%             ener = ener + dist * w_len;
%         end
% 
%         % label orientation bias
% 
%         dx = dx ./ dist;
%         dy = dy ./ dist;
%         if (dx > 0 && dy > 0) 
%             ener = ener + 0 * w_orient; 
%         elseif (dx < 0 && dy > 0) 
%             ener = ener + 1 * w_orient; 
%         elseif (dx < 0 && dy < 0) 
%             ener = ener + 2 * w_orient;
%         else
%             ener = ener + 3 * w_orient;
%         end
% 
%         x21 = lab(index).x;
%         y21 = lab(index).y - lab(index).height + 2.0;
%         x22 = lab(index).x + lab(index).width;
%         y22 = lab(index).y + 2.0;
% 
%         for i = 1:m
%             if (i ~= index)
% 
%                 % penalty for intersection of leader lines
% 
%                 overlap = intersect(anc(index).x, lab(index).x, anc(i).x, lab(i).x,...
%                                     anc(index).y, lab(index).y, anc(i).y, lab(i).y);
%                 if (overlap) 
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
%     end
% end
% 
% 
% function inter = intersect(x1, x2, x3, x4, y1, y2, y3, y4) 
% % returns true if two lines intersect, else false
% % from http:%paulbourke.net/geometry/lineline2d/
% 
%     denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
%     numera = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3);
%     numerb = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3);
% 
%     % Is the intersection along the the segments 
%     mua = numera / denom;
%     mub = numerb / denom;
% 
%     inter = ~(mua < 0 || mua > 1 || mub < 0 || mub > 1);
% 
% end
% 
% function t = cooling_schedule(currT, initialT, nsweeps) 
% % Linear cooling schedule
%     t = (currT - (initialT / nsweeps));
% end
% 



