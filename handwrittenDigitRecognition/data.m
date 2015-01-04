fid = fopen('semion/semeion.data', 'r');
tline = fgetl(fid);
digits = [];
num_digits = 0;
while ischar(tline)
    num_digits = num_digits + 1;
    C = strsplit(tline);
    a = 1;
    b = 16;
    img = [];
    while a <= 241
         var = C(a:b);
         img = [img;var];
         a = a + 16;
         b = b + 16;
    end
    size(img);
    img = str2double(img);
    img = horzcat(zeros(16,3), img);
    img = horzcat(img, zeros(16,3));
    img = vertcat(img, zeros(3,22));
    img = vertcat(zeros(3,22), img);
    if num_digits == 1
        digits = img;
    else
        %digits(:,:,num_digits) = img;
        digits(:,:,num_digits) = img;
    end
    tline = fgetl(fid);
end
disp(num_digits);