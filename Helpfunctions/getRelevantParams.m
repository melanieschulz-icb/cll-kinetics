
function par_positions = getRelevantParams(model)
%     if (1000>model)&& (model > 100)
%         model=model-100;
%     end
    if model == 0 || model==0.5
        %two-compartment model
        par_positions=1:6;
    elseif model == 1
        %two-compartment separation model
        par_positions=[1:6,12];
    elseif model==2
        par_positions=[1,3:6,12];
%     elseif model == 1
%         par_positions=1:6;
%     elseif model == 2
%         par_positions= 1:10;
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

    elseif model == 19
        par_positions=1:12;
    elseif model == 20
        par_positions=[1:9,11,12];
    elseif model == 21
        par_positions=[1:6,8:9,11,12];
    elseif model ==22
        par_positions=[1:6,8:11,12];
    elseif model == 23
        par_positions=[1,3:11,12];
    elseif model ==24
        par_positions=[1,3:6,8:11,12];
    elseif model == 25
        par_positions=[1,3:6,8,9,11,12];
    elseif model == 26
        par_positions=[1,3:9,11,12];
    elseif model == 27
        par_positions=[1,3:5,7:11,12];
    elseif model == 28
        par_positions=[1,3:5,8:11,12];
    elseif model == 29
        par_positions=[1,3:5,7:9,11,12];
    elseif model == 30
        par_positions=[1,3:5,8,9,11,12];
    elseif model == 31
        par_positions=[1:5,7:11,12];
    elseif model == 32
        par_positions=[1:5,8:11,12];
    elseif model == 33
        par_positions=[1:5,7:9,11,12];
    elseif model == 34
        par_positions=[1:5,8,9,11,12];
        
        %no m3
    elseif model == 35
        par_positions=[1:10];
    elseif model == 36
        par_positions=[1:9];
    elseif model == 37
        par_positions=[1:6,8:9];
    elseif model == 38
        par_positions=[1:6,8:10];
    elseif model == 39
        par_positions=[1,3:10];
    elseif model == 40
        par_positions=[1,3:6,8:10];
    elseif model == 41
        par_positions=[1,3:6,8,9];
    elseif model == 42
        par_positions=[1,3:9];
    elseif model == 43
        par_positions=[1,3:5,7:10];
    elseif model == 44
        par_positions=[1,3:5,8:10];
    elseif model == 45
        par_positions=[1,3:5,7:9];
    elseif model == 46
        par_positions=[1,3:5,8,9];
    elseif model == 47
        par_positions=[1:5,7:10];
    elseif model == 48
        par_positions=[1:5,8:10];
    elseif model == 49
        par_positions=[1:5,7:9];
    elseif model == 50
        par_positions=[1:5,8,9];
 
    elseif model == 51
        par_positions=[1:10,12];
    elseif model == 52
        par_positions=[1:9,12];
    elseif model == 53
        par_positions=[1:6,8:9,12];
    elseif model == 54
        par_positions=[1:6,8:10,12];
    elseif model == 55
        par_positions=[1,3:10,12];
    elseif model == 56
        par_positions=[1,3:6,8:10,12];
    elseif model == 57
        par_positions=[1,3:6,8,9,12];
    elseif model == 58
        par_positions=[1,3:9,12];
    elseif model == 59
        par_positions=[1,3:5,7:10,12];
    elseif model == 60
        par_positions=[1,3:5,8:10,12];
    elseif model == 61
        par_positions=[1,3:5,7:9,12];
    elseif model == 62
        par_positions=[1,3:5,8,9,12];
    elseif model == 63
        par_positions=[1:5,7:10,12];
    elseif model == 64
        par_positions=[1:5,8:10,12];
    elseif model == 65
        par_positions=[1:5,7:9,12];
    elseif model == 66
        par_positions=[1:5,8,9,12];
        
        
        %no m1
    elseif model == 67
        par_positions=[1:2,4:11];
    elseif model == 68
        par_positions=[1:2,4:9,11];
    elseif model == 69
        par_positions=[1:2,4:6,8:9,11];
    elseif model == 70
        par_positions=[1:2,4:6,8:11];
    elseif model == 71
        par_positions=[1,4:11];
    elseif model == 72
        par_positions=[1,4:6,8:11];
    elseif model == 73
        par_positions=[1,4:6,8,9,11];
    elseif model == 74
        par_positions=[1,4:9,11];
    elseif model == 75
        par_positions=[1,4:5,7:11];
    elseif model == 76
        par_positions=[1,4:5,8:11];
    elseif model == 77
        par_positions=[1,4:5,7:9,11];
    elseif model == 78
        par_positions=[1,4:5,8,9,11];
    elseif model == 79
        par_positions=[1:2,4:5,7:11];
    elseif model == 80
        par_positions=[1:2,4:5,8:11];
    elseif model == 81
        par_positions=[1:2,4:5,7:9,11];
    elseif model == 82
        par_positions=[1:2,4:5,8,9,11];
     
    elseif model == 83
        par_positions=[1:2,4:11,12];
    elseif model == 84
        par_positions=[1:2,4:9,11,12];
    elseif model == 85
        par_positions=[1:2,4:6,8:9,11,12];
    elseif model == 86
        par_positions=[1:2,4:6,8:11,12];
    elseif model == 87
        par_positions=[1,4:11,12];
    elseif model == 88
        par_positions=[1,4:6,8:11,12];
    elseif model == 89
        par_positions=[1,4:6,8,9,11,12];
    elseif model == 90
        par_positions=[1,4:9,11,12];
    elseif model == 91
        par_positions=[1,4:5,7:11,12];
    elseif model == 92
        par_positions=[1,4:5,8:11,12];
    elseif model == 93
        par_positions=[1,4:5,7:9,11,12];
    elseif model == 94
        par_positions=[1,4:5,8,9,11,12];
    elseif model == 95
        par_positions=[1:2,4:5,7:11,12];
    elseif model == 96
        par_positions=[1:2,4:5,8:11,12];
    elseif model == 97
        par_positions=[1:2,4:5,7:9,11,12];
    elseif model == 98
        par_positions=[1:2,4:5,8,9,11,12];
        
        
    %% no m2    
    elseif model == 99
        par_positions=[1:7,9:11];
    elseif model == 100
        par_positions=[1:7,9,11];
    elseif model == 101
        par_positions=[1,3:7,9:11];
    elseif model == 102
        par_positions=[1,3:7,9,11];
     elseif model == 103
         par_positions=[1,3:5,7,9:11];
    elseif model == 104
        par_positions=[1,3:5,7,9,11];
    elseif model == 105
        par_positions=[1:5,7,9:11];
    elseif model == 106
        par_positions=[1:5,7,9,11];
  
    elseif model == 107
        par_positions=[1:7,9:11,12];
    elseif model == 108
        par_positions=[1:7,9,11,12];
    elseif model == 109
        par_positions=[1,3:7,9:11,12];
    elseif model == 110
        par_positions=[1,3:7,9,11,12];
     elseif model == 111
         par_positions=[1,3:5,7,9:11,12];
    elseif model == 112
        par_positions=[1,3:5,7,9,11,12];
    elseif model == 113
        par_positions=[1:5,7,9:11,12];
    elseif model == 114
        par_positions=[1:5,7,9,11,12];
   
%% only m3

    elseif model == 115
        par_positions=[1,2,4:7,9:11];
    elseif model == 116
        par_positions=[1,2,4:7,9,11];
    elseif model == 117
        par_positions=[1,4:7,9:11];
    elseif model == 118
        par_positions=[1,4:7,9,11];
    elseif model == 119
        par_positions=[1,4:5,7,9:11];
    elseif model == 120
        par_positions=[1,4:5,7,9,11];
    elseif model == 121
        par_positions=[1,2,4:5,7,9:11];
    elseif model == 122
        par_positions=[1,2,4:5,7,9,11];
    elseif model == 123
        par_positions=[1,2,4:7,9:11,12];
    elseif model == 124
        par_positions=[1,2,4:7,9,11,12];
    elseif model == 125
        par_positions=[1,4:7,9:11,12];
    elseif model == 126
        par_positions=[1,4:7,9,11,12];
    elseif model == 127
        par_positions=[1,4:5,7,9:11,12];
    elseif model == 128
        par_positions=[1,4:5,7,9,11,12];
    elseif model == 129
        par_positions=[1,2,4:5,7,9:11,12];
    elseif model == 130
        par_positions=[1,2,4:5,7,9,11,12];

    %% only m
    elseif model == 131
        par_positions=[1:7,9,10];
    elseif model == 132
        par_positions=[1:7,9];
    elseif model == 133
        par_positions=[1:6,9];
    elseif model == 134
        par_positions=[1,3:7,9,10];
    elseif model == 135
        par_positions=[1,3:6,9];
    elseif model == 136
        par_positions=[1,3:7,9];
    elseif model == 137
        par_positions=[1,3:5,7,9,10];
    elseif model == 138
        par_positions=[1,3:5,7,9];
    elseif model == 139
        par_positions=[1,3:5,9];
    elseif model == 140
        par_positions=[1:5,7,9,10];
    elseif model == 141
        par_positions=[1:5,7,9];
    elseif model == 142
        par_positions=[1:5,9];

    elseif model == 143
        par_positions=[1:7,9,10,12];
    elseif model == 144
        par_positions=[1:7,9,12];
    elseif model == 145
        par_positions=[1:6,9,12];
    elseif model == 146
        par_positions=[1,3:7,9,10,12];
    elseif model == 147
        par_positions=[1,3:6,9,12];
    elseif model == 148
        par_positions=[1,3:7,9,12];
    elseif model == 149
        par_positions=[1,3:5,7,9,10,12];
    elseif model == 150
        par_positions=[1,3:5,7,9,12];
    elseif model == 151
        par_positions=[1,3:5,9,12];
    elseif model == 152
        par_positions=[1:5,7,9,10,12];
    elseif model == 153
        par_positions=[1:5,7,9,12];
    elseif model == 154
        par_positions=[1:5,9,12];
        
        
%no m
    elseif model == 155
        par_positions=[1,2,4:7,9,10];
    elseif model == 156
        par_positions=[1,2,4:7,9];
    elseif model == 157
        par_positions=[1,2,4:6,9];
    elseif model == 158
        par_positions=[1,4:5,7,9,10];
    elseif model == 159
        par_positions=[1,4:5,7,9];
    elseif model == 160
        par_positions=[1,4:5,9];
    elseif model == 161
        par_positions=[1,2,4:5,7,9,10];
    elseif model == 162
        par_positions=[1,2,4:5,7,9];
    elseif model == 163
        par_positions=[1,2,4:5,9];
    elseif model == 164
        par_positions=[2,4:7,9,10];
    elseif model == 165
        par_positions=[2,4:7,9];
    elseif model == 166
        par_positions=[2,4:6,9];
    elseif model == 167
        par_positions=[4:5,7,9,10];
    elseif model == 168
        par_positions=[4:5,7,9];
    elseif model == 169
        par_positions=[4:5,9];
    elseif model == 170
        par_positions=[2,4:5,7,9,10];
    elseif model == 171
        par_positions=[2,4:5,7,9];
    elseif model == 172
        par_positions=[2,4:5,9];
    elseif model == 173
        par_positions=[1,2,4:7,9,10,12];
    elseif model == 174
        par_positions=[1,2,4:7,9,12];
    elseif model == 175
        par_positions=[1,2,4:6,9,12];
    elseif model == 176
        par_positions=[1,4:5,7,9,10,12];
    elseif model == 177
        par_positions=[1,4:5,7,9,12];
    elseif model == 178
        par_positions=[1,4:5,9,12];
    elseif model == 179
        par_positions=[1,2,4:5,7,9,10,12];
    elseif model == 180
        par_positions=[1,2,4:5,7,9,12];
    elseif model == 181
        par_positions=[1,2,4:5,9,12];
        
        
%% switch BM and LT
    elseif model == 182
        par_positions=[1     7     8     9     5    10     2     3     4     6    11];
    elseif model == 183
        par_positions=[1     7     8     9     5    10     2     3     4    11];
    elseif model == 184
        par_positions=[1     7     8     9     5    10     3     4    11];
    elseif model == 185
        par_positions=[1     7     8     9     5    10     3     4     6    11];
    elseif model == 186
        par_positions=[1     8     9     5    10     2     3     4     6    11];
    elseif model == 187
        par_positions=[1     8     9     5    10     3     4     6    11];
    elseif model == 188
        par_positions=[1     8     9     5    10     3     4    11];
    elseif model == 189
        par_positions=[1     8     9     5    10     2     3     4    11];
    elseif model == 190
        par_positions=[1     8     9     5     2     3     4     6    11];
    elseif model == 191
        par_positions=[1     8     9     5     3     4     6    11];
    elseif model == 192
        par_positions=[1     8     9     5     2     3     4    11];
    elseif model == 193
        par_positions=[1     8     9     5     3     4    11];
    elseif model == 194
        par_positions=[1     7     8     9     5     2     3     4     6    11];
    elseif model == 195
        par_positions=[1     7     8     9     5     3     4     6    11];
    elseif model == 196
        par_positions=[1     7     8     9     5     2     3     4    11];
    elseif model == 197
        par_positions=[1     7     8     9     5     3     4    11];

    elseif model == 198
        par_positions=[1     7     8     9     5    10     2     3     4     6    11    12];
    elseif model == 199
        par_positions=[1     7     8     9     5    10     2     3     4    11    12];
    elseif model == 200
        par_positions=[1     7     8     9     5    10     3     4    11    12];
    elseif model ==201
        par_positions=[1     7     8     9     5    10     3     4     6    11    12];
    elseif model == 202
        par_positions=[1     8     9     5    10     2     3     4     6    11    12];
    elseif model ==203
        par_positions=[1     8     9     5    10     3     4     6    11    12];
    elseif model == 204
        par_positions=[1     8     9     5    10     3     4    11    12];
    elseif model == 205
        par_positions=[1     8     9     5    10     2     3     4    11    12];
    elseif model == 206
        par_positions=[1     8     9     5     2     3     4     6    11    12];
    elseif model == 207
        par_positions=[1     8     9     5     3     4     6    11    12];
    elseif model == 208
        par_positions=[1     8     9     5     2     3     4    11    12];
    elseif model == 209
        par_positions=[1     8     9     5     3     4    11    12];
    elseif model == 210
        par_positions=[1     7     8     9     5     2     3     4     6    11    12];
    elseif model == 211
        par_positions=[1     7     8     9     5     3     4     6    11    12];
    elseif model == 212
        par_positions=[1     7     8     9     5     2     3     4    11    12];
    elseif model == 213
        par_positions=[1     7     8     9     5     3     4    11    12];
        
        %no m3
    elseif model == 214
        par_positions=[1     7     8     9     5    10     2     3     4     6];
    elseif model == 215
        par_positions=[1     7     8     9     5    10     2     3     4];
    elseif model == 216
        par_positions=[1     7     8     9     5    10     3     4];
    elseif model == 217
        par_positions=[1     7     8     9     5    10     3     4     6];
    elseif model == 218
        par_positions=[1     8     9     5    10     2     3     4     6];
    elseif model == 219
        par_positions=[1     8     9     5    10     3     4     6];
    elseif model == 220
        par_positions=[1     8     9     5    10     3     4];
    elseif model == 221
        par_positions=[1     8     9     5    10     2     3     4];
    elseif model == 222
        par_positions=[1     8     9     5     2     3     4     6];
    elseif model == 223
        par_positions=[1     8     9     5     3     4     6];
    elseif model == 224
        par_positions=[1     8     9     5     2     3     4];
    elseif model == 225
        par_positions=[1     8     9     5     3     4];
    elseif model == 226
        par_positions=[1     7     8     9     5     2     3     4     6];
    elseif model == 227
        par_positions=[1     7     8     9     5     3     4     6];
    elseif model == 228
        par_positions=[1     7     8     9     5     2     3     4];
    elseif model == 229
        par_positions=[1     7     8     9     5     3     4];
 
    elseif model == 230
        par_positions=[1     7     8     9     5    10     2     3     4     6    12];
    elseif model == 231
        par_positions=[1     7     8     9     5    10     2     3     4    12];
    elseif model == 232
        par_positions=[1     7     8     9     5    10     3     4    12];
    elseif model == 233
        par_positions=[1     7     8     9     5    10     3     4     6    12];
    elseif model == 234
        par_positions=[1     8     9     5    10     2     3     4     6    12];
    elseif model == 235
        par_positions=[1     8     9     5    10     3     4     6    12];
    elseif model == 236
        par_positions=[1     8     9     5    10     3     4    12];
    elseif model == 237
        par_positions=[1     8     9     5    10     2     3     4    12];
    elseif model == 238
        par_positions=[1     8     9     5     2     3     4     6    12];
    elseif model == 239
        par_positions=[1     8     9     5     3     4     6    12];
    elseif model == 240
        par_positions=[1     8     9     5     2     3     4    12];
    elseif model == 241
        par_positions=[1     8     9     5     3     4    12];
    elseif model == 242
        par_positions=[1     7     8     9     5     2     3     4     6    12];
    elseif model == 243
        par_positions=[1     7     8     9     5     3     4     6    12];
    elseif model == 244
        par_positions=[1     7     8     9     5     2     3     4    12];
    elseif model == 245
        par_positions=[1     7     8     9     5     3     4    12];
        
        
        %no m1
    elseif model == 246
        par_positions=[1     7     9     5    10     2     3     4     6    11];
    elseif model == 247
        par_positions=[1     7     9     5    10     2     3     4    11];
    elseif model == 248
        par_positions=[1     7     9     5    10     3     4    11];
    elseif model == 249
        par_positions=[1     7     9     5    10     3     4     6    11];
    elseif model == 250
        par_positions=[1     9     5    10     2     3     4     6    11];
    elseif model == 251
        par_positions=[1     9     5    10     3     4     6    11];
    elseif model == 252
        par_positions=[1     9     5    10     3     4    11];
    elseif model == 253
        par_positions=[1     9     5    10     2     3     4    11];
    elseif model == 254
        par_positions=[1     9     5     2     3     4     6    11];
    elseif model == 255
        par_positions=[1     9     5     3     4     6    11];
    elseif model == 256
        par_positions=[1     9     5     2     3     4    11];
    elseif model == 257
        par_positions=[1     9     5     3     4    11];
    elseif model == 258
        par_positions=[1     7     9     5     2     3     4     6    11];
    elseif model == 259
        par_positions=[1     7     9     5     3     4     6    11];
    elseif model == 260
        par_positions=[1     7     9     5     2     3     4    11];
    elseif model == 261
        par_positions=[1     7     9     5     3     4    11];
     
    elseif model == 262
        par_positions=[1     7     9     5    10     2     3     4     6    11    12];
    elseif model == 263
        par_positions=[1     7     9     5    10     2     3     4    11    12];
    elseif model == 264
        par_positions=[1     7     9     5    10     3     4    11    12];
    elseif model == 265
        par_positions=[1     7     9     5    10     3     4     6    11    12];
    elseif model == 266
        par_positions=[1     9     5    10     2     3     4     6    11    12];
    elseif model == 267
        par_positions=[1     9     5    10     3     4     6    11    12];
    elseif model == 268
        par_positions=[1     9     5    10     3     4    11    12];
    elseif model == 269
        par_positions=[1     9     5    10     2     3     4    11    12];
    elseif model == 270
        par_positions=[1     9     5     2     3     4     6    11    12];
    elseif model == 271
        par_positions=[1     9     5     3     4     6    11    12];
    elseif model == 272
        par_positions=[1     9     5     2     3     4    11    12];
    elseif model == 273
        par_positions=[1     9     5     3     4    11    12];
    elseif model == 274
        par_positions=[1     7     9     5     2     3     4     6    11    12];
    elseif model == 275
        par_positions=[1     7     9     5     3     4     6    11    12];
    elseif model == 276
        par_positions=[1     7     9     5     2     3     4    11    12];
    elseif model == 277
        par_positions=[1     7     9     5     3     4    11    12];
        
        
    %% no m2    
    elseif model == 278
        par_positions=[1     7     8     9     5    10     2     4     6    11];
    elseif model == 279
        par_positions=[1     7     8     9     5    10     2     4    11];
    elseif model == 280
        par_positions=[1     8     9     5    10     2     4     6    11];
    elseif model == 281
        par_positions=[1     8     9     5    10     2     4    11];
     elseif model == 282
         par_positions=[1     8     9     5     2     4     6    11];
    elseif model == 283
        par_positions=[1     8     9     5     2     4    11];
    elseif model == 284
        par_positions=[1     7     8     9     5     2     4     6    11];
    elseif model == 285
        par_positions=[1     7     8     9     5     2     4    11];
  
    elseif model == 286
        par_positions=[1     7     8     9     5    10     2     4     6    11    12];
    elseif model == 287
        par_positions=[1     7     8     9     5    10     2     4    11    12];
    elseif model == 288
        par_positions=[1     8     9     5    10     2     4     6    11    12];
    elseif model == 289
        par_positions=[1     8     9     5    10     2     4    11    12];
     elseif model == 290
         par_positions=[1     8     9     5     2     4     6    11    12];
    elseif model == 291
        par_positions=[1     8     9     5     2     4    11    12];
    elseif model == 292
        par_positions=[1     7     8     9     5     2     4     6    11    12];
    elseif model == 293
        par_positions=[1     7     8     9     5     2     4    11    12];
   
%% only m3

    elseif model == 294
        par_positions=[1     7     9     5    10     2     4     6    11];
    elseif model == 295
        par_positions=[1     7     9     5    10     2     4    11];
    elseif model == 296
        par_positions=[1     9     5    10     2     4     6    11];
    elseif model == 297
        par_positions=[1     9     5    10     2     4    11];
    elseif model == 298
        par_positions=[1     9     5     2     4     6    11];
    elseif model == 299
        par_positions=[1     9     5     2     4    11];
    elseif model == 300
        par_positions=[1     7     9     5     2     4     6    11];
    elseif model == 301
        par_positions=[1     7     9     5     2     4    11];
    
    elseif model == 302
        par_positions=[1     7     9     5    10     2     4     6    11    12];
    elseif model == 303
        par_positions=[1     7     9     5    10     2     4    11    12];
    elseif model == 304
        par_positions=[1     9     5    10     2     4     6    11    12];
    elseif model == 305
        par_positions=[1     9     5    10     2     4    11    12];
    elseif model == 306
        par_positions=[1     9     5     2     4     6    11    12];
    elseif model == 307
        par_positions=[1     9     5     2     4    11    12];
    elseif model == 308
        par_positions=[1     7     9     5     2     4     6    11    12];
    elseif model == 309
        par_positions=[1     7     9     5     2     4    11    12];
     
        
 %% switch BM and PB
    elseif model == 310
        par_positions=[7     2    11     4     9     6     1     8     5    12     3];
    elseif model == 311
        par_positions=[7     2    11     4     9     6     1     8     5     3];
    elseif model == 312
        par_positions=[7     2    11     4     9     6     8     5     3];
    elseif model == 313
        par_positions=[7     2    11     4     9     6     8     5    12     3];
    elseif model == 314
        par_positions=[7    11     4     9     6     1     8     5    12     3];
    elseif model == 315
        par_positions=[7    11     4     9     6     8     5    12     3];
    elseif model == 316
        par_positions=[7    11     4     9     6     8     5     3];
    elseif model == 317
        par_positions=[7    11     4     9     6     1     8     5     3];
    elseif model == 318
        par_positions=[7    11     4     9     1     8     5    12     3];
    elseif model == 319
        par_positions=[7    11     4     9     8     5    12     3];
    elseif model == 320
        par_positions=[7    11     4     9     1     8     5     3];
    elseif model == 321
        par_positions=[7    11     4     9     8     5     3];
    elseif model == 322
        par_positions=[7     2    11     4     9     1     8     5    12     3];
    elseif model == 323
        par_positions=[7     2    11     4     9     8     5    12     3];
    elseif model == 324
        par_positions=[7     2    11     4     9     1     8     5     3];
    elseif model == 325
        par_positions=[7     2    11     4     9     8     5     3];

    elseif model == 326
        par_positions=[7     2    11     4     9     6     1     8     5    12     3    10];
    elseif model == 327
        par_positions=[7     2    11     4     9     6     1     8     5     3    10];
    elseif model == 328
        par_positions=[7     2    11     4     9     6     8     5     3    10];
    elseif model ==329
        par_positions=[ 7     2    11     4     9     6     8     5    12     3    10];
    elseif model == 330
        par_positions=[7    11     4     9     6     1     8     5    12     3    10];
    elseif model ==331
        par_positions=[7    11     4     9     6     8     5    12     3    10];
    elseif model == 332
        par_positions=[7    11     4     9     6     8     5     3    10];
    elseif model == 333
        par_positions=[7    11     4     9     6     1     8     5     3    10];
    elseif model == 334
        par_positions=[7    11     4     9     1     8     5    12     3    10];
    elseif model == 335
        par_positions=[7    11     4     9     8     5    12     3    10];
    elseif model == 336
        par_positions=[7    11     4     9     1     8     5     3    10];
    elseif model == 337
        par_positions=[7    11     4     9     8     5     3    10];
    elseif model == 338
        par_positions=[7     2    11     4     9     1     8     5    12     3    10];
    elseif model == 339
        par_positions=[7     2    11     4     9     8     5    12     3    10];
    elseif model == 340
        par_positions=[7     2    11     4     9     1     8     5     3    10];
    elseif model == 341
        par_positions=[7     2    11     4     9     8     5     3    10];
        
        %no m3
    elseif model == 342
        par_positions=[7     2    11     4     9     6     1     8     5    12];
    elseif model == 343
        par_positions=[7     2    11     4     9     6     1     8     5];
    elseif model == 344
        par_positions=[7     2    11     4     9     6     8     5];
    elseif model == 345
        par_positions=[7     2    11     4     9     6     8     5    12];
    elseif model == 346
        par_positions=[7    11     4     9     6     1     8     5    12];
    elseif model == 347
        par_positions=[7    11     4     9     6     8     5    12];
    elseif model == 348
        par_positions=[7    11     4     9     6     8     5];
    elseif model == 349
        par_positions=[7    11     4     9     6     1     8     5];
    elseif model == 350
        par_positions=[7    11     4     9     1     8     5    12];
    elseif model == 351
        par_positions=[7    11     4     9     8     5    12];
    elseif model == 352
        par_positions=[7    11     4     9     1     8     5];
    elseif model == 353
        par_positions=[7    11     4     9     8     5];
    elseif model == 354
        par_positions=[7     2    11     4     9     1     8     5    12];
    elseif model == 355
        par_positions=[7     2    11     4     9     8     5    12];
    elseif model == 356
        par_positions=[7     2    11     4     9     1     8     5];
    elseif model == 357
        par_positions=[7     2    11     4     9     8     5];
 
    elseif model == 358
        par_positions=[7     2    11     4     9     6     1     8     5    12    10];
    elseif model == 359
        par_positions=[7     2    11     4     9     6     1     8     5    10];
    elseif model == 360
        par_positions=[7     2    11     4     9     6     8     5    10];
    elseif model == 361
        par_positions=[7     2    11     4     9     6     8     5    12    10];
    elseif model == 362
        par_positions=[7    11     4     9     6     1     8     5    12    10];
    elseif model == 363
        par_positions=[7    11     4     9     6     8     5    12    10];
    elseif model == 364
        par_positions=[7    11     4     9     6     8     5    10];
    elseif model == 365
        par_positions=[7    11     4     9     6     1     8     5    10];
    elseif model == 366
        par_positions=[7    11     4     9     1     8     5    12    10];
    elseif model == 367
        par_positions=[7    11     4     9     8     5    12    10];
    elseif model == 368
        par_positions=[7    11     4     9     1     8     5    10];
    elseif model == 369
        par_positions=[7    11     4     9     8     5    10];
    elseif model == 370
        par_positions=[7     2    11     4     9     1     8     5    12    10];
    elseif model == 371
        par_positions=[7     2    11     4     9     8     5    12    10];
    elseif model == 372
        par_positions=[7     2    11     4     9     1     8     5    10];
    elseif model == 373
        par_positions=[7     2    11     4     9     8     5    10];
        
        
        %no m1
    elseif model == 374
        par_positions=[7     2     4     9     6     1     8     5    12     3];
    elseif model == 375
        par_positions=[7     2     4     9     6     1     8     5     3];
    elseif model == 376
        par_positions=[7     2     4     9     6     8     5     3];
    elseif model == 377
        par_positions=[7     2     4     9     6     8     5    12     3];
    elseif model == 378
        par_positions=[7     4     9     6     1     8     5    12     3];
    elseif model == 379
        par_positions=[7     4     9     6     8     5    12     3];
    elseif model == 380
        par_positions=[7     4     9     6     8     5     3];
    elseif model == 381
        par_positions=[7     4     9     6     1     8     5     3];
    elseif model == 382
        par_positions=[7     4     9     1     8     5    12     3];
    elseif model == 383
        par_positions=[7     4     9     8     5    12     3];
    elseif model == 384
        par_positions=[7     4     9     1     8     5     3];
    elseif model == 385
        par_positions=[7     4     9     8     5     3];
    elseif model == 386
        par_positions=[7     2     4     9     1     8     5    12     3];
    elseif model == 387
        par_positions=[7     2     4     9     8     5    12     3];
    elseif model == 388
        par_positions=[7     2     4     9     1     8     5     3];
    elseif model == 389
        par_positions=[7     2     4     9     8     5     3];
     
    elseif model == 390
        par_positions=[7     2     4     9     6     1     8     5    12     3    10];
    elseif model == 391
        par_positions=[7     2     4     9     6     1     8     5     3    10];
    elseif model == 392
        par_positions=[7     2     4     9     6     8     5     3    10];
    elseif model == 393
        par_positions=[7     2     4     9     6     8     5    12     3    10];
    elseif model == 394
        par_positions=[7     4     9     6     1     8     5    12     3    10];
    elseif model == 395
        par_positions=[7     4     9     6     8     5    12     3    10];
    elseif model == 396
        par_positions=[7     4     9     6     8     5     3    10];
    elseif model == 397
        par_positions=[7     4     9     6     1     8     5     3    10];
    elseif model == 398
        par_positions=[7     4     9     1     8     5    12     3    10];
    elseif model == 399
        par_positions=[7     4     9     8     5    12     3    10];
    elseif model == 400
        par_positions=[7     4     9     1     8     5     3    10];
    elseif model == 401
        par_positions=[7     4     9     8     5     3    10];
    elseif model == 402
        par_positions=[7     2     4     9     1     8     5    12     3    10];
    elseif model == 403
        par_positions=[7     2     4     9     8     5    12     3    10];
    elseif model == 404
        par_positions=[7     2     4     9     1     8     5     3    10];
    elseif model == 405
        par_positions=[7     2     4     9     8     5     3    10];
        
        
    %% no m2    
    elseif model == 406
        par_positions=[7     2    11     4     9     6     1     5    12     3];
    elseif model == 407
        par_positions=[7     2    11     4     9     6     1     5     3];
    elseif model == 408
        par_positions=[7    11     4     9     6     1     5    12     3];
    elseif model == 409
        par_positions=[7    11     4     9     6     1     5     3];
     elseif model == 410
         par_positions=[7    11     4     9     1     5    12     3];
    elseif model == 411
        par_positions=[7    11     4     9     1     5     3];
    elseif model == 412
        par_positions=[7     2    11     4     9     1     5    12     3];
    elseif model == 413
        par_positions=[7     2    11     4     9     1     5     3];
  
    elseif model == 414
        par_positions=[7     2    11     4     9     6     1     5    12     3    10];
    elseif model == 415
        par_positions=[7     2    11     4     9     6     1     5     3    10];
    elseif model == 416
        par_positions=[7    11     4     9     6     1     5    12     3    10];
    elseif model == 417
        par_positions=[7    11     4     9     6     1     5     3    10];
     elseif model == 418
         par_positions=[7    11     4     9     1     5    12     3    10];
    elseif model == 419
        par_positions=[7    11     4     9     1     5     3    10];
    elseif model == 420
        par_positions=[7     2    11     4     9     1     5    12     3    10];
    elseif model == 421
        par_positions=[7     2    11     4     9     1     5     3    10];
   
%% only m3

    elseif model == 422
        par_positions=[7     2     4     9     6     1     5    12     3];
    elseif model == 423
        par_positions=[7     2     4     9     6     1     5     3];
    elseif model == 424
        par_positions=[7     4     9     6     1     5    12     3];
    elseif model == 425
        par_positions=[7     4     9     6     1     5     3];
    elseif model == 426
        par_positions=[7     4     9     1     5    12     3];
    elseif model == 427
        par_positions=[7     4     9     1     5     3];
    elseif model == 428
        par_positions=[7     2     4     9     1     5    12     3];
    elseif model == 429
        par_positions=[7     2     4     9     1     5     3];
    
    elseif model == 430
        par_positions=[7     2     4     9     6     1     5    12     3    10];
    elseif model == 431
        par_positions=[7     2     4     9     6     1     5     3    10];
    elseif model == 432
        par_positions=[7     4     9     6     1     5    12     3    10];
    elseif model == 433
        par_positions=[7     4     9     6     1     5     3    10];
    elseif model == 434
        par_positions=[7     4     9     1     5    12     3    10];
    elseif model == 435
        par_positions=[7     4     9     1     5     3    10];
    elseif model == 436
        par_positions=[7     2     4     9     1     5    12     3    10];
    elseif model == 437
        par_positions=[7     2     4     9     1     5     3    10];
    
        
%% z->x->y
    elseif model == 438
        par_positions=[2     7    11     9     4    10     1     3     5    12     8];
    elseif model == 439
        par_positions=[2     7    11     9     4    10     1     3     5     8];
    elseif model == 440
        par_positions=[2     7    11     9     4    10     3     5     8];
    elseif model == 441
        par_positions=[2     7    11     9     4    10     3     5    12     8];
    elseif model == 442
        par_positions=[2    11     9     4    10     1     3     5    12     8];
    elseif model == 443
        par_positions=[2    11     9     4    10     3     5    12     8];
    elseif model == 444
        par_positions=[2    11     9     4    10     3     5     8];
    elseif model == 445
        par_positions=[2    11     9     4    10     1     3     5     8];
    elseif model == 446
        par_positions=[2    11     9     4     1     3     5    12     8];
    elseif model == 447
        par_positions=[2    11     9     4     3     5    12     8];
    elseif model == 448
        par_positions=[2    11     9     4     1     3     5     8];
    elseif model == 449
        par_positions=[2    11     9     4     3     5     8];
    elseif model == 450
        par_positions=[2     7    11     9     4     1     3     5    12     8];
    elseif model == 451
        par_positions=[2     7    11     9     4     3     5    12     8];
    elseif model == 452
        par_positions=[2     7    11     9     4     1     3     5     8];
    elseif model == 453
        par_positions=[2     7    11     9     4     3     5     8];

    elseif model == 454
        par_positions=[2     7    11     9     4    10     1     3     5    12     8     6];
    elseif model == 455
        par_positions=[2     7    11     9     4    10     1     3     5     8     6];
    elseif model == 456
        par_positions=[2     7    11     9     4    10     3     5     8     6];
    elseif model ==457
        par_positions=[2     7    11     9     4    10     3     5    12     8     6];
    elseif model == 458
        par_positions=[2    11     9     4    10     1     3     5    12     8     6];
    elseif model ==459
        par_positions=[2    11     9     4    10     3     5    12     8     6];
    elseif model == 460
        par_positions=[2    11     9     4    10     3     5     8     6];
    elseif model == 461
        par_positions=[2    11     9     4    10     1     3     5     8     6];
    elseif model == 462
        par_positions=[2    11     9     4     1     3     5    12     8     6];
    elseif model == 463
        par_positions=[2    11     9     4     3     5    12     8     6];
    elseif model == 464
        par_positions=[2    11     9     4     1     3     5     8     6];
    elseif model == 465
        par_positions=[2    11     9     4     3     5     8     6];
    elseif model == 466
        par_positions=[2     7    11     9     4     1     3     5    12     8     6];
    elseif model == 467
        par_positions=[2     7    11     9     4     3     5    12     8     6];
    elseif model == 468
        par_positions=[2     7    11     9     4     1     3     5     8     6];
    elseif model == 469
        par_positions=[2     7    11     9     4     3     5     8     6];
        
        %no m3
    elseif model == 470
        par_positions=[2     7    11     9     4    10     1     3     5    12];
    elseif model == 471
        par_positions=[2     7    11     9     4    10     1     3     5];
    elseif model == 472
        par_positions=[2     7    11     9     4    10     3     5];
    elseif model == 473
        par_positions=[2     7    11     9     4    10     3     5    12];
    elseif model == 474
        par_positions=[2    11     9     4    10     1     3     5    12];
    elseif model == 475
        par_positions=[2    11     9     4    10     3     5    12];
    elseif model == 476
        par_positions=[2    11     9     4    10     3     5];
    elseif model == 477
        par_positions=[2    11     9     4    10     1     3     5];
    elseif model == 478
        par_positions=[2    11     9     4     1     3     5    12];
    elseif model == 479
        par_positions=[2    11     9     4     3     5    12];
    elseif model == 480
        par_positions=[2    11     9     4     1     3     5];
    elseif model == 481
        par_positions=[2    11     9     4     3     5];
    elseif model == 482
        par_positions=[2     7    11     9     4     1     3     5    12];
    elseif model == 483
        par_positions=[2     7    11     9     4     3     5    12];
    elseif model == 484
        par_positions=[2     7    11     9     4     1     3     5];
    elseif model == 485
        par_positions=[2     7    11     9     4     3     5];
 
    elseif model == 486
        par_positions=[2     7    11     9     4    10     1     3     5    12     6];
    elseif model == 487
        par_positions=[2     7    11     9     4    10     1     3     5     6];
    elseif model == 488
        par_positions=[2     7    11     9     4    10     3     5     6];
    elseif model == 489
        par_positions=[2     7    11     9     4    10     3     5    12     6];
    elseif model == 490
        par_positions=[2    11     9     4    10     1     3     5    12     6];
    elseif model == 491
        par_positions=[2    11     9     4    10     3     5    12     6];
    elseif model == 492
        par_positions=[2    11     9     4    10     3     5     6];
    elseif model == 493
        par_positions=[2    11     9     4    10     1     3     5     6];
    elseif model == 494
        par_positions=[2    11     9     4     1     3     5    12     6];
    elseif model == 495
        par_positions=[2    11     9     4     3     5    12     6];
    elseif model == 496
        par_positions=[2    11     9     4     1     3     5     6];
    elseif model == 497
        par_positions=[2    11     9     4     3     5     6];
    elseif model == 498
        par_positions=[2     7    11     9     4     1     3     5    12     6];
    elseif model == 499
        par_positions=[2     7    11     9     4     3     5    12     6];
    elseif model == 500
        par_positions=[2     7    11     9     4     1     3     5     6];
    elseif model == 501
        par_positions=[2     7    11     9     4     3     5     6];
        
        
        %no m1
    elseif model == 502
        par_positions=[2     7     9     4    10     1     3     5    12     8];
    elseif model == 503
        par_positions=[2     7     9     4    10     1     3     5     8];
    elseif model == 504
        par_positions=[2     7     9     4    10     3     5     8];
    elseif model == 505
        par_positions=[2     7     9     4    10     3     5    12     8];
    elseif model == 506
        par_positions=[2     9     4    10     1     3     5    12     8];
    elseif model == 507
        par_positions=[2     9     4    10     3     5    12     8];
    elseif model == 508
        par_positions=[2     9     4    10     3     5     8];
    elseif model == 509
        par_positions=[2     9     4    10     1     3     5     8];
    elseif model == 510
        par_positions=[2     9     4     1     3     5    12     8];
    elseif model == 511
        par_positions=[2     9     4     3     5    12     8];
    elseif model == 512
        par_positions=[2     9     4     1     3     5     8];
    elseif model == 513
        par_positions=[2     9     4     3     5     8];
    elseif model == 514
        par_positions=[2     7     9     4     1     3     5    12     8];
    elseif model == 515
        par_positions=[2     7     9     4     3     5    12     8];
    elseif model == 516
        par_positions=[2     7     9     4     1     3     5     8];
    elseif model == 517
        par_positions=[2     7     9     4     3     5     8];
     
    elseif model == 518
        par_positions=[2     7     9     4    10     1     3     5    12     8     6];
    elseif model == 519
        par_positions=[2     7     9     4    10     1     3     5     8     6];
    elseif model == 520
        par_positions=[2     7     9     4    10     3     5     8     6];
    elseif model == 521
        par_positions=[2     7     9     4    10     3     5    12     8     6];
    elseif model == 522
        par_positions=[2     9     4    10     1     3     5    12     8     6];
    elseif model == 523
        par_positions=[2     9     4    10     3     5    12     8     6];
    elseif model == 524
        par_positions=[2     9     4    10     3     5     8     6];
    elseif model == 525
        par_positions=[2     9     4    10     1     3     5     8     6];
    elseif model == 526
        par_positions=[2     9     4     1     3     5    12     8     6];
    elseif model == 527
        par_positions=[2     9     4     3     5    12     8     6];
    elseif model == 528
        par_positions=[2     9     4     1     3     5     8     6];
    elseif model == 529
        par_positions=[2     9     4     3     5     8     6];
    elseif model == 530
        par_positions=[2     7     9     4     1     3     5    12     8     6];
    elseif model == 531
        par_positions=[2     7     9     4     3     5    12     8     6];
    elseif model == 532
        par_positions=[2     7     9     4     1     3     5     8     6];
    elseif model == 533
        par_positions=[2     7     9     4     3     5     8     6];
        
        
    %% no m2    
    elseif model == 534
        par_positions=[2     7    11     9     4    10     1     5    12     8];
    elseif model == 535
        par_positions=[2     7    11     9     4    10     1     5     8];
    elseif model == 536
        par_positions=[2    11     9     4    10     1     5    12     8];
    elseif model == 537
        par_positions=[2    11     9     4    10     1     5     8];
     elseif model == 538
         par_positions=[2    11     9     4     1     5    12     8];
    elseif model == 539
        par_positions=[2    11     9     4     1     5     8];
    elseif model == 540
        par_positions=[2     7    11     9     4     1     5    12     8];
    elseif model == 541
        par_positions=[2     7    11     9     4     1     5     8];
  
    elseif model == 542
        par_positions=[2     7    11     9     4    10     1     5    12     8     6];
    elseif model == 543
        par_positions=[2     7    11     9     4    10     1     5     8     6];
    elseif model == 544
        par_positions=[2    11     9     4    10     1     5    12     8     6];
    elseif model == 545
        par_positions=[2    11     9     4    10     1     5     8     6];
     elseif model == 546
         par_positions=[2    11     9     4     1     5    12     8     6];
    elseif model == 547
        par_positions=[2    11     9     4     1     5     8     6];
    elseif model == 548
        par_positions=[2     7    11     9     4     1     5    12     8     6];
    elseif model == 549
        par_positions=[2     7    11     9     4     1     5     8     6];
   
%% only m3

    elseif model == 550
        par_positions=[2     7     9     4    10     1     5    12     8];
    elseif model == 551
        par_positions=[2     7     9     4    10     1     5     8];
    elseif model == 552
        par_positions=[2     9     4    10     1     5    12     8];
    elseif model == 553
        par_positions=[2     9     4    10     1     5     8];
    elseif model == 554
        par_positions=[2     9     4     1     5    12     8];
    elseif model == 555
        par_positions=[2     9     4     1     5     8];
    elseif model == 556
        par_positions=[2     7     9     4     1     5    12     8];
    elseif model == 557
        par_positions=[2     7     9     4     1     5     8];
    
    elseif model == 558
        par_positions=[2     7     9     4    10     1     5    12     8     6];
    elseif model == 559
        par_positions=[2     7     9     4    10     1     5     8     6];
    elseif model == 560
        par_positions=[2     9     4    10     1     5    12     8     6];
    elseif model == 561
        par_positions=[2     9     4    10     1     5     8     6];
    elseif model == 562
        par_positions=[2     9     4     1     5    12     8     6];
    elseif model == 563
        par_positions=[2     9     4     1     5     8     6];
    elseif model == 564
        par_positions=[2     7     9     4     1     5    12     8     6];
    elseif model == 565
        par_positions=[2     7     9     4     1     5     8     6];
        
    
        
    %% y->x->z
    elseif model == 566
        par_positions=[7     1     8     5     9    12     2    11     4     6     3];
    elseif model == 567
        par_positions=[7     1     8     5     9    12     2    11     4     3];
    elseif model == 568
        par_positions=[7     1     8     5     9    12    11     4     3];
    elseif model == 569
        par_positions=[7     1     8     5     9    12    11     4     6     3];
    elseif model == 570
        par_positions=[7     8     5     9    12     2    11     4     6     3];
    elseif model == 571
        par_positions=[7     8     5     9    12    11     4     6     3];
    elseif model == 572
        par_positions=[7     8     5     9    12    11     4     3];
    elseif model == 573
        par_positions=[7     8     5     9    12     2    11     4     3];
    elseif model == 574
        par_positions=[7     8     5     9     2    11     4     6     3];
    elseif model == 575
        par_positions=[7     8     5     9    11     4     6     3];
    elseif model == 576
        par_positions=[7     8     5     9     2    11     4     3];
    elseif model == 577
        par_positions=[7     8     5     9    11     4     3];
    elseif model == 578
        par_positions=[7     1     8     5     9     2    11     4     6     3];
    elseif model == 579
        par_positions=[7     1     8     5     9    11     4     6     3];
    elseif model == 580
        par_positions=[7     1     8     5     9     2    11     4     3];
    elseif model == 581
        par_positions=[7     1     8     5     9    11     4     3];

    elseif model == 582
        par_positions=[7     1     8     5     9    12     2    11     4     6     3    10];
    elseif model == 583
        par_positions=[7     1     8     5     9    12     2    11     4     3    10];
    elseif model == 584
        par_positions=[7     1     8     5     9    12    11     4     3    10];
    elseif model ==585
        par_positions=[7     1     8     5     9    12    11     4     6     3    10];
    elseif model == 586
        par_positions=[7     8     5     9    12     2    11     4     6     3    10];
    elseif model ==587
        par_positions=[7     8     5     9    12    11     4     6     3    10];
    elseif model == 588
        par_positions=[7     8     5     9    12    11     4     3    10];
    elseif model == 589
        par_positions=[7     8     5     9    12     2    11     4     3    10];
    elseif model == 590
        par_positions=[7     8     5     9     2    11     4     6     3    10];
    elseif model == 591
        par_positions=[7     8     5     9    11     4     6     3    10];
    elseif model == 592
        par_positions=[7     8     5     9     2    11     4     3    10];
    elseif model == 593
        par_positions=[7     8     5     9    11     4     3    10];
    elseif model == 594
        par_positions=[7     1     8     5     9     2    11     4     6     3    10];
    elseif model == 595
        par_positions=[7     1     8     5     9    11     4     6     3    10];
    elseif model == 596
        par_positions=[7     1     8     5     9     2    11     4     3    10];
    elseif model == 597
        par_positions=[7     1     8     5     9    11     4     3    10];
        
        %no m3
    elseif model == 598
        par_positions=[7     1     8     5     9    12     2    11     4     6];
    elseif model == 599
        par_positions=[7     1     8     5     9    12     2    11     4];
    elseif model == 600
        par_positions=[7     1     8     5     9    12    11     4];
    elseif model == 601
        par_positions=[7     1     8     5     9    12    11     4     6];
    elseif model == 602
        par_positions=[7     8     5     9    12     2    11     4     6];
    elseif model == 603
        par_positions=[7     8     5     9    12    11     4     6];
    elseif model == 604
        par_positions=[7     8     5     9    12    11     4];
    elseif model == 605
        par_positions=[7     8     5     9    12     2    11     4];
    elseif model == 606
        par_positions=[7     8     5     9     2    11     4     6];
    elseif model == 607
        par_positions=[7     8     5     9    11     4     6];
    elseif model == 608
        par_positions=[7     8     5     9     2    11     4];
    elseif model == 609
        par_positions=[7     8     5     9    11     4];
    elseif model == 610
        par_positions=[7     1     8     5     9     2    11     4     6];
    elseif model == 611
        par_positions=[7     1     8     5     9    11     4     6];
    elseif model == 612
        par_positions=[7     1     8     5     9     2    11     4];
    elseif model == 613
        par_positions=[7     1     8     5     9    11     4];
 
    elseif model == 614
        par_positions=[7     1     8     5     9    12     2    11     4     6    10];
    elseif model == 615
        par_positions=[7     1     8     5     9    12     2    11     4    10];
    elseif model == 616
        par_positions=[7     1     8     5     9    12    11     4    10];
    elseif model == 617
        par_positions=[7     1     8     5     9    12    11     4     6    10];
    elseif model == 618
        par_positions=[7     8     5     9    12     2    11     4     6    10];
    elseif model == 619
        par_positions=[7     8     5     9    12    11     4     6    10];
    elseif model == 620
        par_positions=[7     8     5     9    12    11     4    10];
    elseif model == 621
        par_positions=[7     8     5     9    12     2    11     4    10];
    elseif model == 622
        par_positions=[7     8     5     9     2    11     4     6    10];
    elseif model == 623
        par_positions=[7     8     5     9    11     4     6    10];
    elseif model == 624
        par_positions=[7     8     5     9     2    11     4    10];
    elseif model == 625
        par_positions=[7     8     5     9    11     4    10];
    elseif model == 626
        par_positions=[7     1     8     5     9     2    11     4     6    10];
    elseif model == 627
        par_positions=[7     1     8     5     9    11     4     6    10];
    elseif model == 628
        par_positions=[7     1     8     5     9     2    11     4    10];
    elseif model == 629
        par_positions=[7     1     8     5     9    11     4    10];
        
        
        %no m1
    elseif model == 630
        par_positions=[7     1     5     9    12     2    11     4     6     3];
    elseif model == 631
        par_positions=[7     1     5     9    12     2    11     4     3];
    elseif model == 632
        par_positions=[7     1     5     9    12    11     4     3];
    elseif model == 633
        par_positions=[7     1     5     9    12    11     4     6     3];
    elseif model == 634
        par_positions=[7     5     9    12     2    11     4     6     3];
    elseif model == 635
        par_positions=[7     5     9    12    11     4     6     3];
    elseif model == 636
        par_positions=[7     5     9    12    11     4     3];
    elseif model == 637
        par_positions=[7     5     9    12     2    11     4     3];
    elseif model == 638
        par_positions=[7     5     9     2    11     4     6     3];
    elseif model == 639
        par_positions=[7     5     9    11     4     6     3];
    elseif model == 640
        par_positions=[7     5     9     2    11     4     3];
    elseif model == 641
        par_positions=[7     5     9    11     4     3];
    elseif model == 642
        par_positions=[7     1     5     9     2    11     4     6     3];
    elseif model == 643
        par_positions=[7     1     5     9    11     4     6     3];
    elseif model == 644
        par_positions=[7     1     5     9     2    11     4     3];
    elseif model == 645
        par_positions=[7     1     5     9    11     4     3];
     
    elseif model == 646
        par_positions=[7     1     5     9    12     2    11     4     6     3    10];
    elseif model == 647
        par_positions=[7     1     5     9    12     2    11     4     3    10];
    elseif model == 648
        par_positions=[7     1     5     9    12    11     4     3    10];
    elseif model == 649
        par_positions=[7     1     5     9    12    11     4     6     3    10];
    elseif model == 650
        par_positions=[7     5     9    12     2    11     4     6     3    10];
    elseif model == 651
        par_positions=[7     5     9    12    11     4     6     3    10];
    elseif model == 652
        par_positions=[7     5     9    12    11     4     3    10];
    elseif model == 653
        par_positions=[7     5     9    12     2    11     4     3    10];
    elseif model == 654
        par_positions=[7     5     9     2    11     4     6     3    10];
    elseif model == 655
        par_positions=[7     5     9    11     4     6     3    10];
    elseif model == 656
        par_positions=[7     5     9     2    11     4     3    10];
    elseif model == 657
        par_positions=[7     5     9    11     4     3    10];
    elseif model == 658
        par_positions=[7     1     5     9     2    11     4     6     3    10];
    elseif model == 659
        par_positions=[7     1     5     9    11     4     6     3    10];
    elseif model == 660
        par_positions=[7     1     5     9     2    11     4     3    10];
    elseif model == 661
        par_positions=[7     1     5     9    11     4     3    10];
        
        
    %% no m2    
    elseif model == 662
        par_positions=[7     1     8     5     9    12     2     4     6     3];
    elseif model == 663
        par_positions=[7     1     8     5     9    12     2     4     3];
    elseif model == 664
        par_positions=[7     8     5     9    12     2     4     6     3];
    elseif model == 665
        par_positions=[7     8     5     9    12     2     4     3];
     elseif model == 666
         par_positions=[7     8     5     9     2     4     6     3];
    elseif model == 667
        par_positions=[7     8     5     9     2     4     3];
    elseif model == 668
        par_positions=[7     1     8     5     9     2     4     6     3];
    elseif model == 669
        par_positions=[7     1     8     5     9     2     4     3];
  
    elseif model == 670
        par_positions=[7     1     8     5     9    12     2     4     6     3    10];
    elseif model == 671
        par_positions=[7     1     8     5     9    12     2     4     3    10];
    elseif model == 672
        par_positions=[7     8     5     9    12     2     4     6     3    10];
    elseif model == 673
        par_positions=[7     8     5     9    12     2     4     3    10];
     elseif model == 674
         par_positions=[7     8     5     9     2     4     6     3    10];
    elseif model == 675
        par_positions=[7     8     5     9     2     4     3    10];
    elseif model == 676
        par_positions=[7     1     8     5     9     2     4     6     3    10];
    elseif model == 677
        par_positions=[7     1     8     5     9     2     4     3    10];
   
%% only m3

    elseif model == 678
        par_positions=[7     1     5     9    12     2     4     6     3];
    elseif model == 679
        par_positions=[7     1     5     9    12     2     4     3];
    elseif model == 680
        par_positions=[7     5     9    12     2     4     6     3];
    elseif model == 681
        par_positions=[7     5     9    12     2     4     3];
    elseif model == 682
        par_positions=[7     5     9     2     4     6     3];
    elseif model == 683
        par_positions=[7     5     9     2     4     3];
    elseif model == 684
        par_positions=[7     1     5     9     2     4     6     3];
    elseif model == 685
        par_positions=[7     1     5     9     2     4     3];
    
    elseif model == 686
        par_positions=[7     1     5     9    12     2     4     6     3    10];
    elseif model == 687
        par_positions=[7     1     5     9    12     2     4     3    10];
    elseif model == 688
        par_positions=[7     5     9    12     2     4     6     3    10];
    elseif model == 689
        par_positions=[7     5     9    12     2     4     3    10];
    elseif model == 690
        par_positions=[7     5     9     2     4     6     3    10];
    elseif model == 691
        par_positions=[7     5     9     2     4     3    10];
    elseif model == 692
        par_positions=[7     1     5     9     2     4     6     3    10];
    elseif model == 693
        par_positions=[7     1     5     9     2     4     3    10];
    
    
        
    %% switch PB and LT
    elseif model == 694
        par_positions=[2     1     3     5     4    12     7    11     9    10     8];
    elseif model == 695
        par_positions=[2     1     3     5     4    12     7    11     9     8];
    elseif model == 696
        par_positions=[2     1     3     5     4    12    11     9     8];
    elseif model == 697
        par_positions=[2     1     3     5     4    12    11     9    10     8];
    elseif model == 698
        par_positions=[2     3     5     4    12     7    11     9    10     8];
    elseif model == 699
        par_positions=[2     3     5     4    12    11     9    10     8];
    elseif model == 700
        par_positions=[2     3     5     4    12    11     9     8];
    elseif model == 701
        par_positions=[2     3     5     4    12     7    11     9     8];
    elseif model == 702
        par_positions=[2     3     5     4     7    11     9    10     8];
    elseif model == 703
        par_positions=[2     3     5     4    11     9    10     8];
    elseif model == 704
        par_positions=[2     3     5     4     7    11     9     8];
    elseif model == 705
        par_positions=[2     3     5     4    11     9     8];
    elseif model == 706
        par_positions=[2     1     3     5     4     7    11     9    10     8];
    elseif model == 707
        par_positions=[2     1     3     5     4    11     9    10     8];
    elseif model == 708
        par_positions=[2     1     3     5     4     7    11     9     8];
    elseif model == 709
        par_positions=[2     1     3     5     4    11     9     8];

    elseif model == 710
        par_positions=[2     1     3     5     4    12     7    11     9    10     8     6];
    elseif model == 711
        par_positions=[2     1     3     5     4    12     7    11     9     8     6];
    elseif model == 712
        par_positions=[2     1     3     5     4    12    11     9     8     6];
    elseif model == 713
        par_positions=[2     1     3     5     4    12    11     9    10     8     6];
    elseif model == 714
        par_positions=[2     3     5     4    12     7    11     9    10     8     6];
    elseif model ==715
        par_positions=[2     3     5     4    12    11     9    10     8     6];
    elseif model == 716
        par_positions=[2     3     5     4    12    11     9     8     6];
    elseif model == 717
        par_positions=[2     3     5     4    12     7    11     9     8     6];
    elseif model == 718
        par_positions=[2     3     5     4     7    11     9    10     8     6];
    elseif model == 719
        par_positions=[ 2     3     5     4    11     9    10     8     6];
    elseif model == 720
        par_positions=[2     3     5     4     7    11     9     8     6];
    elseif model == 721
        par_positions=[2     3     5     4    11     9     8     6];
    elseif model == 722
        par_positions=[2     1     3     5     4     7    11     9    10     8     6];
    elseif model == 723
        par_positions=[2     1     3     5     4    11     9    10     8     6];
    elseif model == 724
        par_positions=[2     1     3     5     4     7    11     9     8     6];
    elseif model == 725
        par_positions=[2     1     3     5     4    11     9     8     6];
        
        %no m3
    elseif model == 726
        par_positions=[2     1     3     5     4    12     7    11     9    10];
    elseif model == 727
        par_positions=[2     1     3     5     4    12     7    11     9];
    elseif model == 728
        par_positions=[2     1     3     5     4    12    11     9];
    elseif model == 729
        par_positions=[2     1     3     5     4    12    11     9    10];
    elseif model == 730
        par_positions=[2     3     5     4    12     7    11     9    10];
    elseif model == 731
        par_positions=[2     3     5     4    12    11     9    10];
    elseif model == 732
        par_positions=[2     3     5     4    12    11     9];
    elseif model == 733
        par_positions=[2     3     5     4    12     7    11     9];
    elseif model == 734
        par_positions=[2     3     5     4     7    11     9    10];
    elseif model == 735
        par_positions=[2     3     5     4    11     9    10];
    elseif model == 736
        par_positions=[2     3     5     4     7    11     9];
    elseif model == 737
        par_positions=[2     3     5     4    11     9];
    elseif model == 738
        par_positions=[2     1     3     5     4     7    11     9    10];
    elseif model == 739
        par_positions=[2     1     3     5     4    11     9    10];
    elseif model == 740
        par_positions=[2     1     3     5     4     7    11     9];
    elseif model == 741
        par_positions=[2     1     3     5     4    11     9];
 
    elseif model == 742
        par_positions=[2     1     3     5     4    12     7    11     9    10     6];
    elseif model == 743
        par_positions=[2     1     3     5     4    12     7    11     9     6];
    elseif model == 744
        par_positions=[2     1     3     5     4    12    11     9     6];
    elseif model == 745
        par_positions=[2     1     3     5     4    12    11     9    10     6];
    elseif model == 746
        par_positions=[2     3     5     4    12     7    11     9    10     6];
    elseif model == 747
        par_positions=[2     3     5     4    12    11     9    10     6];
    elseif model == 748
        par_positions=[2     3     5     4    12    11     9     6];
    elseif model == 749
        par_positions=[2     3     5     4    12     7    11     9     6];
    elseif model == 750
        par_positions=[2     3     5     4     7    11     9    10     6];
    elseif model == 751
        par_positions=[2     3     5     4    11     9    10     6];
    elseif model == 752
        par_positions=[2     3     5     4     7    11     9     6];
    elseif model == 753
        par_positions=[2     3     5     4    11     9     6];
    elseif model == 754
        par_positions=[2     1     3     5     4     7    11     9    10     6];
    elseif model == 755
        par_positions=[2     1     3     5     4    11     9    10     6];
    elseif model == 756
        par_positions=[2     1     3     5     4     7    11     9     6];
    elseif model == 757
        par_positions=[2     1     3     5     4    11     9     6];
        
        
        %no m1
    elseif model == 758
        par_positions=[2     1     5     4    12     7    11     9    10     8];
    elseif model == 759
        par_positions=[2     1     5     4    12     7    11     9     8];
    elseif model == 760
        par_positions=[2     1     5     4    12    11     9     8];
    elseif model == 761
        par_positions=[2     1     5     4    12    11     9    10     8];
    elseif model == 762
        par_positions=[2     5     4    12     7    11     9    10     8];
    elseif model == 763
        par_positions=[2     5     4    12    11     9    10     8];
    elseif model == 764
        par_positions=[2     5     4    12    11     9     8];
    elseif model == 765
        par_positions=[2     5     4    12     7    11     9     8];
    elseif model == 766
        par_positions=[2     5     4     7    11     9    10     8];
    elseif model == 767
        par_positions=[2     5     4    11     9    10     8];
    elseif model == 768
        par_positions=[2     5     4     7    11     9     8];
    elseif model == 769
        par_positions=[2     5     4    11     9     8];
    elseif model == 770
        par_positions=[2     1     5     4     7    11     9    10     8];
    elseif model == 771
        par_positions=[2     1     5     4    11     9    10     8];
    elseif model == 772
        par_positions=[2     1     5     4     7    11     9     8];
    elseif model == 773
        par_positions=[2     1     5     4    11     9     8];
     
    elseif model == 774
        par_positions=[2     1     5     4    12     7    11     9    10     8     6];
    elseif model == 775
        par_positions=[2     1     5     4    12     7    11     9     8     6];
    elseif model == 776
        par_positions=[2     1     5     4    12    11     9     8     6];
    elseif model == 777
        par_positions=[2     1     5     4    12    11     9    10     8     6];
    elseif model == 778
        par_positions=[2     5     4    12     7    11     9    10     8     6];
    elseif model == 779
        par_positions=[2     5     4    12    11     9    10     8     6];
    elseif model == 780
        par_positions=[2     5     4    12    11     9     8     6];
    elseif model == 781
        par_positions=[2     5     4    12     7    11     9     8     6];
    elseif model == 782
        par_positions=[2     5     4     7    11     9    10     8     6];
    elseif model == 783
        par_positions=[2     5     4    11     9    10     8     6];
    elseif model == 784
        par_positions=[2     5     4     7    11     9     8     6];
    elseif model == 785
        par_positions=[2     5     4    11     9     8     6];
    elseif model == 786
        par_positions=[2     1     5     4     7    11     9    10     8     6];
    elseif model == 787
        par_positions=[2     1     5     4    11     9    10     8     6];
    elseif model == 788
        par_positions=[2     1     5     4     7    11     9     8     6];
    elseif model == 789
        par_positions=[2     1     5     4    11     9     8     6];
        
        
    %% no m2    
    elseif model == 790
        par_positions=[2     1     3     5     4    12     7     9    10     8];
    elseif model == 791
        par_positions=[2     1     3     5     4    12     7     9     8];
    elseif model == 792
        par_positions=[2     3     5     4    12     7     9    10     8];
    elseif model == 793
        par_positions=[2     3     5     4    12     7     9     8];
     elseif model == 794
         par_positions=[2     3     5     4     7     9    10     8];
    elseif model == 795
        par_positions=[2     3     5     4     7     9     8];
    elseif model == 796
        par_positions=[2     1     3     5     4     7     9    10     8];
    elseif model == 797
        par_positions=[2     1     3     5     4     7     9     8];
  
    elseif model == 798
        par_positions=[2     1     3     5     4    12     7     9    10     8     6];
    elseif model == 799
        par_positions=[2     1     3     5     4    12     7     9     8     6];
    elseif model == 800
        par_positions=[2     3     5     4    12     7     9    10     8     6];
    elseif model == 801
        par_positions=[2     3     5     4    12     7     9     8     6];
     elseif model == 802
         par_positions=[ 2     3     5     4     7     9    10     8     6];
    elseif model == 803
        par_positions=[2     3     5     4     7     9     8     6];
    elseif model == 804
        par_positions=[2     1     3     5     4     7     9    10     8     6];
    elseif model == 805
        par_positions=[2     1     3     5     4     7     9     8     6];
   
%% only m3

    elseif model == 806
        par_positions=[2     1     5     4    12     7     9    10     8];
    elseif model == 807
        par_positions=[2     1     5     4    12     7     9     8];
    elseif model == 808
        par_positions=[2     5     4    12     7     9    10     8];
    elseif model == 809
        par_positions=[2     5     4    12     7     9     8];
    elseif model == 810
        par_positions=[2     5     4     7     9    10     8];
    elseif model == 811
        par_positions=[2     5     4     7     9     8];
    elseif model == 812
        par_positions=[2     1     5     4     7     9    10     8];
    elseif model == 813
        par_positions=[2     1     5     4     7     9     8];
    
    elseif model == 814
        par_positions=[2     1     5     4    12     7     9    10     8     6];
    elseif model == 815
        par_positions=[2     1     5     4    12     7     9     8     6];
    elseif model == 816
        par_positions=[2     5     4    12     7     9    10     8     6];
    elseif model == 817
        par_positions=[2     5     4    12     7     9     8     6];
    elseif model == 818
        par_positions=[2     5     4     7     9    10     8     6];
    elseif model == 819
        par_positions=[2     5     4     7     9     8     6];
    elseif model == 820
        par_positions=[2     1     5     4     7     9    10     8     6];
    elseif model == 821
        par_positions=[2     1     5     4     7     9     8     6];    
        
        
 %% adapted model, with having constant compartments     
 
    elseif model == 1003
        par_positions=1:12;
    elseif model == 1004    
        par_positions =[1     2     3     4     5     6     8     9    10    11    12];
    elseif model == 1005
        par_positions =[1     3     4     5     6     7     8     9    10    11    12];
    elseif model == 1006
        par_positions =[1     3     4     5     6     8     9    10    11    12];
    elseif model == 1007
        par_positions =[1     2     3     4     5     6     7     8     9    10    12];
    elseif model == 1008
        par_positions =[1     2     3     4     5     6     8     9    10    12];
    elseif model == 1009
        par_positions =[1     3     4     5     6     7     8     9    10    12];
    elseif model == 1010
        par_positions =[1     3     4     5     6     8     9    10    12];
    elseif model == 1011
        par_positions =[1     2     4     5     6     7     8     9    10    11    12];
    elseif model == 1012
        par_positions =[1     2     4     5     6     8     9    10    11    12];
    elseif model == 1013
        par_positions =[1     4     5     6     7     8     9    10    11    12];
    elseif model == 1014
        par_positions =[1     4     5     6     8     9    10    11    12];
    elseif model == 1015
        par_positions =[1     2     3     4     5     6     7     9    10    11    12];
    elseif model == 1016
        par_positions =[1     3     4     5     6     7     9    10    11    12];
    elseif model == 1017
        par_positions =[1     2     4     5     6     7     9    10    11    12];
    elseif model == 1018
        par_positions =[1     4     5     6     7     9    10    11    12];
    elseif model == 1019
        par_positions =[1     2     3     4     5     6     7     9    10    12];
    elseif model == 1020
        par_positions =[1     3     4     5     6     7     9    10    12];
    elseif model == 1021
        par_positions =[1     2     4     5     6     7     9    10    12];
        
%% switch BM and LT
    elseif model == 1022    
        par_positions =[1     7     8     9     5    10     2     3     4     6    11    12];
    elseif model == 1023
        par_positions =[1     7     8     9     5    10     3     4     6    11    12];
    elseif model == 1024
        par_positions =[1     8     9     5    10     2     3     4     6    11    12];
    elseif model == 1025
        par_positions =[1     8     9     5    10     3     4     6    11    12];
    elseif model == 1026
        par_positions =[1     7     8     9     5    10     2     3     4     6    12];
    elseif model == 1027
        par_positions =[1     7     8     9     5    10     3     4     6    12];
    elseif model == 1028
        par_positions =[1     8     9     5    10     2     3     4     6    12];
    elseif model == 1029
        par_positions =[1     8     9     5    10     3     4     6    12];
    elseif model == 1030
        par_positions =[1     7     9     5    10     2     3     4     6    11    12];
    elseif model == 1031
        par_positions =[1     7     9     5    10     3     4     6    11    12];
    elseif model == 1032
        par_positions =[1     9     5    10     2     3     4     6    11    12];
    elseif model == 1033
        par_positions =[1     9     5    10     3     4     6    11    12];
    elseif model == 1034
        par_positions =[1     7     8     9     5    10     2     4     6    11    12];
    elseif model == 1035
        par_positions =[ 1     8     9     5    10     2     4     6    11    12];
    elseif model == 1036
        par_positions =[1     7     9     5    10     2     4     6    11    12];
    elseif model == 1037
        par_positions =[1     9     5    10     2     4     6    11    12];
        
%% switch BM and PB
    elseif model == 1038
        par_positions =[7     2    11     4     9     6     1     8     5    12     3    10];
    elseif model == 1039
        par_positions =[7     2    11     4     9     6     8     5    12     3    10];
    elseif model == 1040
        par_positions =[7    11     4     9     6     1     8     5    12     3    10];
    elseif model == 1041
        par_positions =[7    11     4     9     6     8     5    12     3    10];
    elseif model == 1042
        par_positions =[7     2    11     4     9     6     1     8     5    12    10];
    elseif model == 1043
        par_positions =[7     2    11     4     9     6     8     5    12    10];
    elseif model == 1044
        par_positions =[7    11     4     9     6     1     8     5    12    10];
    elseif model == 1045
        par_positions =[7    11     4     9     6     8     5    12    10];
    elseif model == 1046
        par_positions =[7     2     4     9     6     1     8     5    12     3    10];
    elseif model == 1047
        par_positions =[ 7     2     4     9     6     8     5    12     3    10];
    elseif model == 1048
        par_positions =[7     4     9     6     1     8     5    12     3    10];
    elseif model == 1049
        par_positions =[7     4     9     6     8     5    12     3    10];
    elseif model == 1050
        par_positions =[7     2    11     4     9     6     1     5    12     3    10];
    elseif model == 1051
        par_positions =[7    11     4     9     6     1     5    12     3    10];
    elseif model == 1052
        par_positions =[7     2     4     9     6     1     5    12     3    10];
    elseif model == 1053
        par_positions =[7     4     9     6     1     5    12     3    10];

%% z->x->y
    elseif model == 1054
        par_positions =[2     7    11     9     4    10     1     3     5    12     8     6];
    elseif model == 1055
        par_positions =[2     7    11     9     4    10     3     5    12     8     6];
    elseif model == 1056
        par_positions =[2    11     9     4    10     1     3     5    12     8     6];
    elseif model == 1057
        par_positions =[2    11     9     4    10     3     5    12     8     6];
    elseif model == 1058
        par_positions =[2     7    11     9     4    10     1     3     5    12     6];
    elseif model == 1059
        par_positions =[2     7    11     9     4    10     3     5    12     6];
    elseif model == 1060
        par_positions =[2    11     9     4    10     1     3     5    12     6];
    elseif model == 1061
        par_positions =[2    11     9     4    10     3     5    12     6];
    elseif model == 1062
        par_positions =[2     7     9     4    10     1     3     5    12     8     6];
    elseif model == 1063
        par_positions =[2     7     9     4    10     3     5    12     8     6];
    elseif model == 1064
        par_positions =[2     9     4    10     1     3     5    12     8     6];
    elseif model == 1065
        par_positions =[2     9     4    10     3     5    12     8     6];
    elseif model == 1066
        par_positions =[2     7    11     9     4    10     1     5    12     8     6];
    elseif model == 1067
        par_positions =[2    11     9     4    10     1     5    12     8     6];
    elseif model == 1068
        par_positions =[2     7     9     4    10     1     5    12     8     6];
    elseif model == 1069
        par_positions =[2     9     4    10     1     5    12     8     6];

%% y->x->z
    elseif model == 1070
        par_positions =[7     1     8     5     9    12     2    11     4     6     3    10];
    elseif model == 1071
        par_positions =[7     1     8     5     9    12    11     4     6     3    10];
    elseif model == 1072
        par_positions =[7     8     5     9    12     2    11     4     6     3    10];
    elseif model == 1073
        par_positions =[7     8     5     9    12    11     4     6     3    10];
    elseif model == 1074
        par_positions =[7     1     8     5     9    12     2    11     4     6    10];
    elseif model == 1075
        par_positions =[7     1     8     5     9    12    11     4     6    10];
    elseif model == 1076
        par_positions =[7     8     5     9    12     2    11     4     6    10];
    elseif model == 1077
        par_positions =[7     8     5     9    12    11     4     6    10];
    elseif model == 1078
        par_positions =[7     1     5     9    12     2    11     4     6     3    10];
    elseif model == 1079
        par_positions =[7     1     5     9    12    11     4     6     3    10];
    elseif model == 1080
        par_positions =[7     5     9    12     2    11     4     6     3    10];
    elseif model == 1081
        par_positions =[7     5     9    12    11     4     6     3    10];
    elseif model == 1082
        par_positions =[7     1     8     5     9    12     2     4     6     3    10];
    elseif model == 1083
        par_positions =[ 7     8     5     9    12     2     4     6     3    10];
    elseif model == 1084
        par_positions =[7     1     5     9    12     2     4     6     3    10];
    elseif model == 1085
        par_positions =[7     5     9    12     2     4     6     3    10];

%% switch PB and LT
    elseif model == 1086
        par_positions =[2     1     3     5     4    12     7    11     9    10     8     6];
    elseif model == 1087
        par_positions =[2     1     3     5     4    12    11     9    10     8     6];
    elseif model == 1088
        par_positions =[2     3     5     4    12     7    11     9    10     8     6];
    elseif model == 1089
        par_positions =[2     3     5     4    12    11     9    10     8     6];
    elseif model == 1090
        par_positions =[2     1     3     5     4    12     7    11     9    10     6];
    elseif model == 1091
        par_positions =[2     1     3     5     4    12    11     9    10     6];
    elseif model == 1092
        par_positions =[2     3     5     4    12     7    11     9    10     6];
    elseif model == 1093
        par_positions =[2     3     5     4    12    11     9    10     6];
    elseif model == 1094
        par_positions =[ 2     1     5     4    12     7    11     9    10     8     6];
    elseif model == 1095
        par_positions =[2     1     5     4    12    11     9    10     8     6];
    elseif model == 1096
        par_positions =[2     5     4    12     7    11     9    10     8     6];
    elseif model == 1097
        par_positions =[ 2     5     4    12    11     9    10     8     6];
    elseif model == 1098
        par_positions =[2     1     3     5     4    12     7     9    10     8     6];
    elseif model == 1099
        par_positions =[2     3     5     4    12     7     9    10     8     6];
    elseif model == 1100
        par_positions =[2     1     5     4    12     7     9    10     8     6];
    elseif model == 1101
        par_positions =[2     5     4    12     7     9    10     8     6 ];       
    elseif 2003 <= model <=2101  
        par_with_zh=getRelevantParams(model-1000);
        names=getVarNames(model-1000);

        ind_zh=find(names(getRelevantParams(model-1000))=="BM_healthy");
%        ind_zh=10;
        par_positions=par_with_zh([1:ind_zh-1,ind_zh+1:end]);
        
    else
        error("model not specified");
    end        
end