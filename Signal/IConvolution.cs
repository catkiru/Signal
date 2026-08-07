namespace MatrixConvolution;

public interface IConvolution
{
    bool DoConvolution(int[,] originalMatrix, int[,] kernel, out int[,] resultMatrix);
}