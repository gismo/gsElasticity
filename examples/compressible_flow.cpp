#include <gsCore/gsTemplateTools.h>
#include <gsElasticity/gsCompressibleFlow.h>
#include <gsIO/gsWriteParaview.h>
#include <gsNurbs/gsNurbsCreator.h>
#include <gsSolver/gsIterativeSolver.h>
#include <iostream>

using namespace gismo;

int main()
{
    // 定义维度
    constexpr int dim = 2;

    // 创建一个2D单位正方形的几何对象 (例如：单位正方形)
    gsMultiPatch<> geometry;
    gsNurbsCreator<>::createUnitSquare(geometry);

    // 设置边界条件对象
    gsBoundaryConditions<> bcInfo;
    
    // 设置流体速度和压力的多分片对象
    gsMultiPatch<> velocity, pressure;

    // 初始化为常数
    gsFunctionExpr<> initVelocity("0.1", "0.1", 2);
    gsFunctionExpr<> initPressure("1.0", 2);

    velocity.addPatch(initVelocity);
    pressure.addPatch(initPressure);

    // 创建流体装配器对象
    gsCompressibleFlowAssembler<> assembler(velocity, pressure, bcInfo, geometry, false);

    // 设置系统参数
    gsMatrix<> solVector;
    assembler.assemble(solVector, {});

    // 调用 assemble 装配系统
    assembler.assemble(true);

    // 调用正性检查函数
    assembler.checkPositivity();

    // 将结果写入Paraview文件以进行可视化
    gsWriteParaview(geometry, "compressible_flow_result");

    std::cout << "Compressible flow simulation complete!" << std::endl;
    return 0;
}