namespace MatrixConvolution;

public class ConvolutionChangeSize(string inputPath) : IConvolution
{
    public bool DoConvolution(int[,] originalMatrix, int[,] kernel, out int[,] resultMatrix)
    {
        int originalWidth = originalMatrix.GetLength(1);
        int i = 0;
        int[,] prev = originalMatrix;
        int[,] current = null;
        resultMatrix = null;

        while (prev.GetLength(1) > originalWidth / 3)
        {
            int ah = prev.GetLength(0);
            int aw = prev.GetLength(1);
            current = new int[ah, aw];
            Array.Copy(prev, current, prev.Length);
            var next = new int[ah - 2, aw - 2];

            for (int y = 1; y < ah - 1; y++)
            {
                for (int x = 1; x < aw - 1; x++)
                {
                    var sum =
                        prev[y - 1, x - 1] * kernel[0, 0] +
                        prev[y - 1, x - 0] * kernel[0, 1] +
                        prev[y - 1, x + 1] * kernel[0, 2] +
                        prev[y - 0, x - 1] * kernel[1, 0] +
                        prev[y - 0, x - 0] * kernel[1, 1] +
                        prev[y - 0, x + 1] * kernel[1, 2] +
                        prev[y + 1, x - 1] * kernel[2, 0] +
                        prev[y + 1, x - 0] * kernel[2, 1] +
                        prev[y + 1, x + 1] * kernel[2, 2];
                    current[y, x] = sum > 0 ? 1 : 0;
                    next[y - 1, x - 1] = current[y, x];
                }
            }

            Tools.SaveMatrixAsBmp(current, inputPath + $".sizeless.{i++}.bmp");
            if (Tools.CompareMatrix(current, prev))
            {
                resultMatrix = current;
                return true;
            }

            prev = next;
        }

        return false;
    }
}