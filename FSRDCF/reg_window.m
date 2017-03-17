function [ W reg_window_s] = reg_window(w_sz, target_sz )

 reg_scale = 0.8* target_sz;
 unitary = sqrt(prod(w_sz));

wrg = -(w_sz(1))/2:(w_sz(1)-2)/2;
wcg = -(w_sz(2))/2:(w_sz(2)-2)/2;

[wrs, wcs] = ndgrid(wrg, wcg);

reg_window = (3 - 0.1) * (abs(wrs/reg_scale(1)).^2 + abs(wcs/reg_scale(2)).^2) + 0.1;
    
% compute the DFT and enforce sparsity
reg_window_dft = fft2(reg_window) / prod(w_sz);%
reg_window_dft_sep = cat(3, real(reg_window_dft), imag(reg_window_dft));
reg_window_dft_sep(abs(reg_window_dft_sep) < 0.05 * max(abs(reg_window_dft_sep(:)))) = 0;
reg_window_dft = reg_window_dft_sep(:,:,1) + 1i*reg_window_dft_sep(:,:,2);
    
% do the inverse transform, correct window minimum
reg_window_sparse = real(ifft2(reg_window_dft));
reg_window_dft(1,1) = reg_window_dft(1,1) -  prod(w_sz)* min(reg_window_sparse(:)) + 0.1;%prod(w_sz)
reg_window = real(ifft2(reg_window_dft));
% construct the regukarization window
regf = ifft2(reg_window)*unitary;%
regf = real(regf);
%
regf(abs(regf) < 10^(-10)) = 0;
%  [o, ~] = sort(abs(regf(:)),'descend');
%  regf(abs(regf) < 0.05*max(abs(regf(:)))) = 0;
% reg_window_s = real(fft2(real(regf)))/unitary;
% m = min((reg_window_s(:)));
% regf(1,1) = regf(1,1) - m*unitary+ 0.0001*unitary;%
reg_window_s = real(fft2(real(regf)))/unitary;%
% surf(real(reg_window_s));
w = sparse(real(reg_window_dft));%real(reg/unitaryf
% w = w.^2;
% W = gallery('circul',w(:));
% % temp = eye(size(w,1)-1);
% % temp = [temp zeros(size(temp,1),1)];
% % shiftR = [[zeros(1,size(temp,1)) 1];temp];
% % temp = eye(size(w,2)-1);
% % temp = [zeros(size(temp,2),1) temp];
% % shiftC = [temp; [1 zeros(1,size(temp,2)-1)]];
% % tempC = w;
% % C0 =w;
% % k=1;
% % for i = 1:size(w,2)
% %     if(i>1)
% %         C0 = tempC*shiftC;
% %         tempC = C0;
% %     end
% %     R = C0;
% %     for j =1:size(w,1)
% %         if(j>1)
% %             R = shiftR*C0;
% %             C0 = R;
% %         end
% %         W(k,:) = R(:)';
% %         k = k +1;
% %     end
% % end
% % 
W = cconvmtx2(w)';
W = W'*W;
end

