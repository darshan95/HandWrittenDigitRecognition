function [imgs labels] = readsemeion(imgFile, readDigits, offset)
    fid = fopen(imgFile, 'r');
    tline = fgetl(fid);
    line_num = 1;
    while line_num ~= offset
        tline = fgetl(fid);
        line_num = line_num + 1;
    end
    digits = [];
    num_digits = 0;
    label = [];
    count_digits = 0;
    while count_digits <= readDigits
        count_digits = count_digits + 1;
        num_digits = num_digits + 1;
        disp(tline)
        class(tline)
        %tline = num2str(tline);
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
        %for label
        temp_label_arr = C(257:266);
        temp_label_arr = str2double(temp_label_arr);
        index = 1;
        index_val = 0;
        while temp_label_arr(index) ~= 1
            index = index + 1;
        end
        curr_label = index - 1;
        label = [label;curr_label];
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
    imgs = digits;
    labels = label;
end
