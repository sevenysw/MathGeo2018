function[P] = patch_decenter(Q,options)
or_size = size(Q);
options.null = 0;
w = getoptions(options,'w',sqrt(size(Q,1)));
m = getoptions(options,'m',size(Q,2));
center_pos = getoptions(options,'center_pos',ceil(w^2/2));
linear = getoptions(options,'linear',0);

if not(linear)
    cx = options.orig(1,:);
    cy = options.orig(2,:);
    if (size(Q,3) ==1)
        Q = reshape(Q, [w w m]);
    end
    P = zeros(w,w,m);
    for k = 1:m
        patch = Q(:,:,k);
        %         [~,max_pos] = max(patch(:));
        %         [x,y] = ind2sub([w w],max_pos);
        dx = (center_pos(1,k) - cx);
        dy = (center_pos(2,k) - cy);
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
        
        patch2 = patch(:,ty(end:-1:1));
        patch2 = patch2(tx(end:-1:1),:);
        patch2(:) = patch2(end:-1:1);
        P(:,:,k) = patch2;
    end
    P = reshape(P,or_size);
else
    if size(Q,3)>1
        Q = reshape(Q,[w^2,m]);
    end
    P = zeros(size(Q));
    for k = 1:m
        inv = 0;
        patch = Q(:,k);
        center = centers(k);
        dpos = center_pos - center;
        if dpos>0
            patch = patch(end:-1:1);
            inv = 1;
            dpos = -dpos;
        end
        if dpos<0
            patch2 = [patch(1+-dpos:end);patch(1:-dpos)];
            if inv
                patch2 = patch2(end:-1:1);
            end
        else
            patch2 = patch;
        end
        P(:,k) = patch2;
    end
    P = reshape(P,or_size);
end