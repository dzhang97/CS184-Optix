from scipy import misc
import numpy as np

def shift_image(image, x, y):
    temp = np.zeros_like(image, dtype=np.float64)
    if x < 0:
        temp[:,:x] = image[:,-x:]
    elif x > 0:
        temp[:,x:] = image[:,:-x]
    else:
        temp = image
    ret = np.zeros_like(temp, dtype=np.float64)
    if y < 0:
        ret[:y] = temp[-y:]
    elif y > 0:
        ret[y:] = temp[:-y]
    else:
        ret = temp
    return ret;

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

    def refocus(self, f):
        image = np.zeros(self.images[0].shape)
        shifty = -(f*(self.n-1)/2)
        for j in range(self.n):
            shiftx = -(f*(self.n-1)/2)
            for i in range(self.n):
                image += shift_image(self.images[j*self.n + i], shiftx, shifty)
                shiftx += f
            shifty += f
        return np.round(image / self.n**2).astype(np.uint8)

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    lf = LightField("images/example_image_", 3, 480, 360)
    # plt.figure()
    # plt.imshow(lf.get_refocused(0.5))
    plt.figure()
<<<<<<< HEAD
    temp = lf.refocus(5)
    plt.imshow(temp)
=======
    plt.imshow(lf.get_refocused(0))
    plt.figure()
    plt.imshow(lf.get_refocused(-0.5))
    plt.figure()
    plt.imshow(lf.get_refocused(-0.5))
    plt.figure()
    plt.imshow(lf.get_refocused(-1.0))
    plt.figure()
    plt.imshow(lf.get_refocused(-1.5))
>>>>>>> 5396e67d59cdbe609473ec8a4ce9f01fecadc1e6
    plt.show()