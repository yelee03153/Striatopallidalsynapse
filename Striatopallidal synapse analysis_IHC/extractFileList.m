function image_list = extractFileList(folder_name)

list = dir(folder_name); % ?΄λ―Έμ? ???₯? ?΄??? ??Ό λ¦¬μ€?Έ μΆμΆ

image_list = [];

for ii = 3:length(list) % 1, 2λ²μ§Έ κ°μ? ?? ?? λΆ?λΆ?(?­?) κ·Έλ? 3λΆ??° ?΄λ―Έμ? ??Ό?. κ°μ? ?΄?? ?????Ό κ°μ λΉΌμΌ??κΉ? -1?¨
    image_list = [image_list; string(list(ii).name)];
end

end