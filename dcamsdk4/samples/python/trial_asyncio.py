#%%
import asyncio, numpy as np
from skimage import io
event = asyncio.Event()
%matplotlib qt
async def get_user_input():
  global user_input
  user_input = await asyncio.get_event_loop().run_in_executor(None, input, "Enter something (or press Enter to skip): ")
  event.set()  # Signal user input is available

async def print_numbers():
  while True:
    for i in range(1, 11):
      frame = np.random.randn(100,100)
      io.imshow(frame)

      # Check for user input after each iteration
      user_input = await get_user_input()
      if user_input:
        return  # Exit the coroutine if input is received

async def main():
  # Create tasks for printing and user input
  print_task = asyncio.create_task(print_numbers())
  input_task = asyncio.create_task(get_user_input())

  # Wait for both tasks to complete or one to be cancelled
  await asyncio.gather(*[print_task, input_task])

  # Print any received user input
  if user_input:
    print(f"User entered: {user_input}")

if __name__ == "__main__":
  await main()
