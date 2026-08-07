namespace MatrixConvolution;

public class ConvolutionFixedSize(string inputPath) : IConvolution
{
    public bool DoConvolution(int[,] originalMatrix, int[,] kernel, out int[,] resultMatrix)
    {
        int originalWidth = originalMatrix.GetLength(1);
        int i = 0;
        int[,] prev = originalMatrix;
        int[,] current = null;
        resultMatrix = null;

        while (i < originalWidth / 3)
        {
            int ah = prev.GetLength(0);
            int aw = prev.GetLength(1);
            current = new int[ah, aw];
            Array.Copy(prev, current, prev.Length);

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
                }
            }

            Tools.SaveMatrixAsBmp(current, inputPath + $".fixed.{i}.bmp");
            i++;
         
            if (true)
            {
                resultMatrix = current;
                return true;
            }

            prev = current;
        }

        return false;
    }
}