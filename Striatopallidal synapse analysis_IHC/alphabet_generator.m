function output_matrix = alphabet_generator()
    output_matrix = [];
    for ii = 1:10
        if ii>1
            str = string(char(64+ii-1));
        else
            str = "";
        end
        for jj = 1:26
            output_matrix = [output_matrix str+string(char(64+jj))];
        end

    end
end