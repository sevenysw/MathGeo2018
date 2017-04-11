function[Q,centers] = patch_center(P,options)

options.null = 0;
w = getoptions(options,'w',sqrt(size(P,1)));
m = getoptions(options,'m',size(P,2));
or_size = size(P);
linear = getoptions(options,'linear',0);
if not(linear)
    centers = zeros(2,m);
    if size(P,3) ==1
        P = reshape(P, [w w m]);
    end
    Q = zeros(size(P));
    for k = 1:m
        patch = P(:,:,k);
        center_pos = options.center_pos(:,k);
        [~,max_pos] = max(patch(:));
        [x,y] = ind2sub([w w],max_pos);
        centers(:,k) = [x,y];
        dx = center_pos - x;
        dy = center_pos - y;
        if dx==1
            x0 = w;
        else
            x0 = mod(1-dx,w);
        end
        if dy ==1
            y0 = w;
        else
            y0 = mod(1-dy,w);
        end
        tx = [x0:w,1:x0-1];
        ty = [y0:w,1:y0-1];
        
        patch = patch(tx,:);
        patch = patch(:,ty);
        Q(:,:,k) = patch;
    end
    Q = reshape(Q,[or_size]);
else
    centers = zeros(1,m);
    if size(P,2)>1
        P = reshape(P,[w*w,m]);
    end
    Q = zeros(size(P));
    for k = 1:m
        center_pos = options.center_pos(k);
        inv = 0;
        patch = P(:,k);
        [~,max_pos] = max(patch);
        centers(k) = max_pos;
        dpos = center_pos - max_pos;
        if dpos>0
            patch = patch(end:-1:1);
            inv = 1;
            dpos = -dpos;
        end
        if dpos<0
            patch2 = [patch(1-dpos:end);patch(1:-dpos)];
            if inv
                patch2 = patch2(end:-1:1);
            end
        else
            patch2 = patch;
        end
        Q(:,k) = patch2;
    end
    Q = reshape(Q,or_size);
end