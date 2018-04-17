from scipy import misc
import numpy as np

class LightField:
    
    def __init__(self, path_prefix, n, width=480, height=360):
        self.images = [misc.imread(path_prefix + str(i) + ".png") for i in range(n**2)]
        self.width = width
        self.height = height
        self.n = n
        self.pos = {}
        middle = np.array([int(self.n/2), int(self.n/2)])
        for i in range(n**2):
            self.pos[i] = np.array([i % self.n, i / self.n]) - middle

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

    def get_value(self, i, x, y):
        if (x >= self.width or y >= self.height):
            return np.array([0,0,0])
        return self.images[i][int(y)][int(x)]

    def sample(self, x, y, f):
        color = np.array([0,0,0])
        for i in range(self.n ** 2):
            pos = self.pos[i]
            x_i = x + pos[0]*f
            y_i = y + pos[1]*f
            color += self.get_value(i, x_i, y_i)
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
    lf = LightField("images/example_image_", 3, 945, 705)
    # plt.figure()
    # plt.imshow(lf.get_refocused(0.5))
    plt.figure()
    plt.imshow(lf.get_refocused(-4))
    plt.show()