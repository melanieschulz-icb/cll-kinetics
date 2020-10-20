function par_positions = getRelevantParams(model)
    if model == 1
        par_positions=1:6;
    elseif model == 2
        par_positions= 1:10;
    elseif model == 3
        par_positions=1:11;
    elseif model == 4
        par_positions=[1:9,11];
    elseif model == 5
        par_positions=[1:6,8:9,11];
    elseif model == 6
        par_positions=[1:6,8:11];
    elseif model == 7
        par_positions=[1,3:11];
    elseif model == 8
        par_positions=[1,3:6,8:11];
    elseif model == 9
        par_positions=[1,3:6,8,9,11];
    elseif model == 10
        par_positions=[1,3:9,11];
    elseif model == 11
        par_positions=[1,3:5,7:11];
    elseif model == 12
        par_positions=[1,3:5,8:11];
    elseif model == 13
        par_positions=[1,3:5,7:9,11];
    elseif model == 14
        par_positions=[1,3:5,8,9,11];
    elseif model == 15
        par_positions=[1:5,7:11];
    elseif model == 16
        par_positions=[1:5,8:11];
    elseif model == 17
        par_positions=[1:5,7:9,11];
    elseif model == 18
        par_positions=[1:5,8,9,11];
        %no m3
    elseif model == 19
        par_positions=[1:10];
    elseif model == 20
        par_positions=[1:9];
    elseif model == 21
        par_positions=[1:6,8:9];
    elseif model == 22
        par_positions=[1:6,8:10];
    elseif model == 23
        par_positions=[1,3:10];
    elseif model == 24
        par_positions=[1,3:6,8:10];
    elseif model == 25
        par_positions=[1,3:6,8,9];
    elseif model == 26
        par_positions=[1,3:9];
    elseif model == 27
        par_positions=[1,3:5,7:10];
    elseif model == 28
        par_positions=[1,3:5,8:10];
    elseif model == 29
        par_positions=[1,3:5,7:9];
    elseif model == 30
        par_positions=[1,3:5,8,9];
    elseif model == 31
        par_positions=[1:5,7:10];
    elseif model == 32
        par_positions=[1:5,8:10];
    elseif model == 33
        par_positions=[1:5,7:9];
    elseif model == 34
        par_positions=[1:5,8,9];
        %no m1
    elseif model == 35
        par_positions=[1:2,4:11];
    elseif model == 36
        par_positions=[1:2,4:9,11];
    elseif model == 37
        par_positions=[1:2,4:6,8:9,11];
    elseif model == 38
        par_positions=[1:2,4:6,8:11];
    elseif model == 39
        par_positions=[1,4:11];
    elseif model == 40
        par_positions=[1,4:6,8:11];
    elseif model == 41
        par_positions=[1,4:6,8,9,11];
    elseif model == 42
        par_positions=[1,4:9,11];
    elseif model == 43
        par_positions=[1,4:5,7:11];
    elseif model == 44
        par_positions=[1,4:5,8:11];
    elseif model == 45
        par_positions=[1,4:5,7:9,11];
    elseif model == 46
        par_positions=[1,4:5,8,9,11];
    elseif model == 47
        par_positions=[1:2,4:5,7:11];
    elseif model == 48
        par_positions=[1:2,4:5,8:11];
    elseif model == 49
        par_positions=[1:2,4:5,7:9,11];
    elseif model == 50
        par_positions=[1:2,4:5,8,9,11];
    %no m2    
    elseif model == 51
        par_positions=[1:7,9:11];
    elseif model == 52
        par_positions=[1:7,9,11];
    elseif model == 53
        par_positions=[1:6,9,11];
    elseif model == 54
        par_positions=[1:6,9:11];
    elseif model == 55
        par_positions=[1,3:7,9:11];
    elseif model == 56
        par_positions=[1,3:6,9:11];
    elseif model == 57
        par_positions=[1,3:6,9,11];
    elseif model == 58
        par_positions=[1,3:7,9,11];
    elseif model == 59
        par_positions=[1,3:5,7,9:11];
    elseif model == 60
        par_positions=[1,3:5,9:11];
    elseif model == 61
        par_positions=[1,3:5,7,9,11];
    elseif model == 62
        par_positions=[1,3:5,9,11];
    elseif model == 63
        par_positions=[1:5,7,9:11];
    elseif model == 64
        par_positions=[1:5,9:11];
    elseif model == 65
        par_positions=[1:5,7,9,11];
    elseif model == 66
        par_positions=[1:5,9,11];
    else
        par_positions=[];
    end
end