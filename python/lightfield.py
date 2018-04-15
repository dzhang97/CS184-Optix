from scipy import misc
import numpy as np

class LightField:
    
    def __init__(self, path_prefix, n, width=480, height=360):
        self.images = [misc.imread(path_prefix + str(i) + ".png") for i in range(n**2)]
        self.width = width
        self.height = height
        self.n = n

    def get_square(self):
        square = []
        for h in range(self.n * self.height):
            row = []
            for w in range(self.n * self.width):
                i = (w % self.n) + self.n * (h % self.n)
                row.append(self.get_value(i, w/self.n, h/self.n))
            square.append(row)
        return square
        
    def get_side_by_side(self):
        side_by_side = []
        for h in range(self.height):
            row = []
            for w in range(self.width):
                for i in range(self.n**2):
                    row.append(self.images[i][h][w])
            side_by_side.append(row)
        return side_by_side

    def get_pos(self, i):
        middle = np.array([self.n/2 - 1, self.n/2 - 1])
        pos = np.array([i % self.n, i / self.n]) - middle
        return pos
        # if pos[1] > pos[0]:
        #     return -np.linalg.norm(pos)
        # elif pos[1] < pos[0]:
        #     return np.linalg.norm(pos)
        # else
        #     return 0

    def get_value(self, i, x, y):
        if (x >= self.width or y >= self.height):
            return np.array([0,0,0])
        return self.images[i][int(y)][int(x)]

    def sample(self, x, y, f):
        color = np.array([0,0,0])
        for i in range(self.n ** 2):
            pos = self.get_pos(i)
            x += pos[0]*f
            y += pos[1]*f
            color += self.get_value(i, x, y)
        return color / (self.n ** 2)

    def get_refocused(self, f):
        image = []
        for h in range(self.height):
            row = []
            for w in range(self.width):
                row.append(self.sample(w, h, f))
            image.append(row)
        return image



if __name__ == '__main__':
    from matplotlib import pyplot as plt
    lf = LightField("images/dragon", 3, 480, 360)
    plt.imshow(lf.get_refocused(1))
    plt.show()