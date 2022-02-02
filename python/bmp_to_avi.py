import cv2
import numpy as np
import os

# Parameters to control
image_folder = 'cropped_dimer3_shrink'
remove_kw = 'cropped_dimer3_'
fps = 10
video_name = image_folder + '.avi'
# ==================
frame_dictionary = {}
# Rename all files in specified directory
for img in os.listdir(image_folder):
    temp = img.replace(remove_kw, "")
    if img != temp:
        os.rename(image_folder + "/" + img, image_folder + "/" + temp)
    img = temp
    frame_dictionary[img.replace('.bmp', '')] = img


images = []
# Permute the images to the correct order
for i in range(len(frame_dictionary)):
    images.append(frame_dictionary[str(i+1)])

# Sort image list

frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape
    
video = cv2.VideoWriter(video_name, 0, fps, (width, height))

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))
    
cv2.destroyAllWindows()
video.release()