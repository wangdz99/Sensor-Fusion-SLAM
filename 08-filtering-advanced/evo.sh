# run evo evaluation:
# a. fused 没有输入运动模型  输出评估结果，并以zip的格式存储:
evo_ape kitti ground_truth.txt fused.txt -r full --plot --plot_mode xy  --save_results ./fused.zip
# b. fused_vel 速度观测  输出评估结果，并以zip的格式存储:
evo_ape kitti ground_truth.txt fused.txt -r full --plot --plot_mode xy  --save_results ./fused_vel.zip
# c. fused_cons 运动约束伪观测  输出评估结果，并以zip的格式存储:
evo_ape kitti ground_truth.txt fused.txt -r full --plot --plot_mode xy  --save_results ./fused_cons.zip
#e. 比较 laser  fused  一并比较评估
evo_res  *.zip --use_filenames -p    
